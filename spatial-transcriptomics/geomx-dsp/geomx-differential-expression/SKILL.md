---
name: bio-spatial-transcriptomics-geomx-differential-expression
description: Differential expression analysis for GeoMx DSP data using linear mixed-effect models (LMM). Account for tissue subsampling and compare regions within slides or disease states between slides. Use when comparing gene expression between GeoMx segments.
tool_type: r
primary_tool: GeomxTools
---

# GeoMx Differential Expression

Differential expression using linear mixed-effect models for GeoMx data.

## Required Libraries

```r
library(GeomxTools)
library(ggplot2)
library(ggrepel)
library(reshape2)
```

## Linear Mixed Model Overview

GeoMx uses LMM to account for tissue subsampling (multiple ROIs per slide):

| Comparison Type | Model | Random Effect |
|----------------|-------|---------------|
| **Within-slide** (e.g., glomerulus vs tubule) | `~ region + (1 + region \| slide)` | Random slope + intercept |
| **Between-slide** (e.g., DKD vs healthy) | `~ disease + (1 \| slide)` | Random intercept only |

## Within-Slide Comparison (Random Slope Model)

For comparing features that co-exist within the same tissue (e.g., different anatomical regions).

```r
# Prepare data
pData(target_demoData)$testRegion <-
    factor(pData(target_demoData)$region, c("glomerulus", "tubule"))
pData(target_demoData)[["slide"]] <-
    factor(pData(target_demoData)[["slide name"]])

# Create log-transformed expression
assayDataElement(target_demoData, elt = "log_q") <-
    assayDataApply(target_demoData, 2, FUN = log, base = 2, elt = "q_norm")

# Run LMM for each disease status separately
results <- c()
for(status in c("DKD", "normal")) {
    ind <- pData(target_demoData)$class == status

    # Formula: ~ testRegion (fixed) + (1 + testRegion | slide) (random)
    # Random slope accounts for region effects varying by slide
    mixedOutmc <- mixedModelDE(
        target_demoData[, ind],
        elt = "log_q",
        modelFormula = ~ testRegion + (1 + testRegion | slide),
        groupVar = "testRegion",
        nCores = parallel::detectCores(),
        multiCore = FALSE
    )

    # Format results
    r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
    tests <- rownames(r_test)
    r_test <- as.data.frame(r_test)
    r_test$Contrast <- tests
    r_test$Gene <- unlist(lapply(colnames(mixedOutmc),
                                 rep, nrow(mixedOutmc["lsmeans", ][[1]])))
    r_test$Subset <- status
    r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
    r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate",
                         "Pr(>|t|)", "FDR")]
    results <- rbind(results, r_test)
}

# View results
head(results)
#   Gene  Subset  Contrast               Estimate  Pr(>|t|)  FDR
#   ACTA2 DKD     glomerulus - tubule    -1.066    0.0065    0.112
#   PDHA1 normal  glomerulus - tubule    -1.542    0.0000    0.000
```

## Between-Slide Comparison (Random Intercept Model)

For comparing disease states or treatments that are tissue-specific.

```r
# Prepare data
pData(target_demoData)$testClass <-
    factor(pData(target_demoData)$class, c("normal", "DKD"))

# Run LMM for each region separately
results2 <- c()
for(region in c("glomerulus", "tubule")) {
    ind <- pData(target_demoData)$region == region

    # Formula: ~ testClass (fixed) + (1 | slide) (random)
    # Random intercept accounts for slide-to-slide variation
    mixedOutmc <- mixedModelDE(
        target_demoData[, ind],
        elt = "log_q",
        modelFormula = ~ testClass + (1 | slide),
        groupVar = "testClass",
        nCores = parallel::detectCores(),
        multiCore = FALSE
    )

    # Format results
    r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
    tests <- rownames(r_test)
    r_test <- as.data.frame(r_test)
    r_test$Contrast <- tests
    r_test$Gene <- unlist(lapply(colnames(mixedOutmc),
                                 rep, nrow(mixedOutmc["lsmeans", ][[1]])))
    r_test$Subset <- region
    r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
    r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate",
                         "Pr(>|t|)", "FDR")]
    results2 <- rbind(results2, r_test)
}
```

## Interpreting Results

```r
# View top genes
top_genes <- subset(results, FDR < 0.05 & abs(Estimate) > 1)
top_genes <- top_genes[order(top_genes$FDR), ]
head(top_genes, 20)

# Contrast interpretation:
# "glomerulus - tubule" with Estimate = 2.5
#   -> Gene is 2^2.5 = 5.7-fold higher in glomerulus
# "glomerulus - tubule" with Estimate = -1.5
#   -> Gene is 2^1.5 = 2.8-fold higher in tubule
```

## Volcano Plot

```r
# Categorize by significance
results$Color <- "NS or FC < 0.5"
results$Color[results$`Pr(>|t|)` < 0.05] <- "P < 0.05"
results$Color[results$FDR < 0.05] <- "FDR < 0.05"
results$Color[results$FDR < 0.001] <- "FDR < 0.001"
results$Color[abs(results$Estimate) < 0.5] <- "NS or FC < 0.5"
results$Color <- factor(results$Color,
                        levels = c("NS or FC < 0.5", "P < 0.05",
                                   "FDR < 0.05", "FDR < 0.001"))

# Select top genes for labeling
results$invert_P <- (-log10(results$`Pr(>|t|)`)) * sign(results$Estimate)
top_g <- c()
for(cond in unique(results$Subset)) {
    ind <- results$Subset == cond
    top_g <- c(top_g,
               results[ind, 'Gene'][order(results[ind, 'invert_P'],
                                         decreasing = TRUE)[1:15]],
               results[ind, 'Gene'][order(results[ind, 'invert_P'],
                                         decreasing = FALSE)[1:15]])
}
top_g <- unique(top_g)

# Plot volcano
ggplot(results,
       aes(x = Estimate, y = -log10(`Pr(>|t|)`),
           color = Color, label = Gene)) +
    geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
    geom_hline(yintercept = -log10(0.05), lty = "dashed") +
    geom_point() +
    labs(x = "log2(Fold Change)", y = "-log10(P-value)") +
    scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                  `FDR < 0.05` = "lightblue",
                                  `P < 0.05` = "orange2",
                                  `NS or FC < 0.5` = "gray")) +
    geom_text_repel(data = subset(results, Gene %in% top_g & FDR < 0.001),
                    size = 3, max.overlaps = 50) +
    theme_bw() +
    facet_wrap(~Subset, scales = "free_y")
```

## MA Plot

```r
# Add mean expression
results$MeanExp <- rowMeans(assayDataElement(target_demoData, elt = "q_norm"))

# Plot
ggplot(results, aes(x = MeanExp, y = Estimate,
                    size = -log10(`Pr(>|t|)`),
                    color = Color, label = Gene)) +
    geom_hline(yintercept = c(0.5, -0.5), lty = "dashed") +
    scale_x_continuous(trans = "log2") +
    geom_point(alpha = 0.5) +
    labs(y = "log2(Fold Change)", x = "Mean Expression") +
    scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                  `FDR < 0.05` = "lightblue",
                                  `P < 0.05` = "orange2",
                                  `NS or FC < 0.5` = "gray")) +
    theme_bw() +
    facet_wrap(~Subset)
```

## Plot Individual Genes

```r
# Violin plot for specific gene
gene_of_interest <- "NPHS1"

ggplot(pData(target_demoData),
       aes(x = region, fill = region,
           y = as.numeric(assayDataElement(target_demoData[gene_of_interest, ],
                                elt = "q_norm")))) +
    geom_violin() +
    geom_jitter(width = .2) +
    labs(y = paste(gene_of_interest, "Expression")) +
    scale_y_continuous(trans = "log2") +
    facet_wrap(~class) +
    theme_bw()
```

## Heatmap of Significant Genes

```r
library(pheatmap)

# Select genes with FDR < 0.001
GOI <- unique(subset(results, FDR < 0.001)$Gene)

# Plot heatmap
pheatmap(log2(assayDataElement(target_demoData[GOI, ], elt = "q_norm")),
         scale = "row",
         show_rownames = FALSE,
         show_colnames = FALSE,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         cutree_cols = 2,
         cutree_rows = 2,
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = pData(target_demoData)[, c("region", "class")])
```

## Export Results

```r
# Write results to CSV
write.csv(results, file = "geomx_de_results.csv", row.names = FALSE)

# Write significant genes only
sig_genes <- subset(results, FDR < 0.05 & abs(Estimate) > 1)
write.csv(sig_genes, file = "geomx_de_significant.csv", row.names = FALSE)
```

## Complete DE Workflow

```r
library(GeomxTools)
library(ggplot2)

# Load normalized data
target_demoData <- readRDS("geomx_normalized.rds")

# Prepare log-transformed data
assayDataElement(target_demoData, elt = "log_q") <-
    assayDataApply(target_demoData, 2, FUN = log, base = 2, elt = "q_norm")

# Set up factors
pData(target_demoData)$testRegion <-
    factor(pData(target_demoData)$region, c("glomerulus", "tubule"))
pData(target_demoData)$slide <-
    factor(pData(target_demoData)$`slide name`)

# Run within-slide DE (glomerulus vs tubule)
results <- c()
for(status in c("DKD", "normal")) {
    ind <- pData(target_demoData)$class == status
    mixedOutmc <- mixedModelDE(
        target_demoData[, ind],
        elt = "log_q",
        modelFormula = ~ testRegion + (1 + testRegion | slide),
        groupVar = "testRegion",
        nCores = parallel::detectCores(),
        multiCore = FALSE
    )

    r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
    r_test <- as.data.frame(r_test)
    r_test$Contrast <- rownames(r_test)
    r_test$Gene <- unlist(lapply(colnames(mixedOutmc),
                                 rep, nrow(mixedOutmc["lsmeans", ][[1]])))
    r_test$Subset <- status
    r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
    results <- rbind(results, r_test[, c("Gene", "Subset", "Contrast",
                                         "Estimate", "Pr(>|t|)", "FDR")])
}

# Export
write.csv(results, file = "de_results.csv", row.names = FALSE)

# Summarize
sig_summary <- table(
    Subset = results$Subset,
    Significant = results$FDR < 0.05
)
print(sig_summary)
```

## Related Skills

- geomx-normalization - Normalize before DE
- differential-expression/deseq2-basics - Alternative DE approach
- data-visualization/volcano-customization - Volcano plot styling
