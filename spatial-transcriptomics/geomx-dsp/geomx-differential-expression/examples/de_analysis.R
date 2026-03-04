#!/usr/bin/env Rscript
# Differential expression analysis for GeoMx data

library(GeomxTools)
library(ggplot2)

# Load normalized data
# target_demoData <- readRDS("geomx_normalized.rds")

cat("Differential Expression Analysis\n")
cat("Data:", nrow(target_demoData), "genes x", ncol(target_demoData), "segments\n\n")

# Prepare log-transformed data
cat("Step 1: Creating log2-transformed expression...\n")
assayDataElement(target_demoData, elt = "log_q") <-
    assayDataApply(target_demoData, 2, FUN = log, base = 2, elt = "q_norm")

# Set up factors
cat("Step 2: Preparing factors for analysis...\n")
pData(target_demoData)$testRegion <-
    factor(pData(target_demoData)$region, c("glomerulus", "tubule"))
pData(target_demoData)$slide <-
    factor(pData(target_demoData)$`slide name`)

# Run within-slide comparison: glomerulus vs tubule
cat("\nStep 3: Running LMM for glomerulus vs tubule...\n")
cat("  Using random slope model: ~ region + (1 + region | slide)\n\n")

results <- c()
for(status in c("DKD", "normal")) {
    cat("  Analyzing", status, "samples...\n")
    ind <- pData(target_demoData)$class == status

    # Run mixed model
    # modelFormula: random slope (1 + testRegion | slide) accounts for region effects varying by slide
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
    r_test <- as.data.frame(r_test)
    r_test$Contrast <- rownames(r_test)
    r_test$Gene <- unlist(lapply(colnames(mixedOutmc),
                                 rep, nrow(mixedOutmc["lsmeans", ][[1]])))
    r_test$Subset <- status
    r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
    r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate",
                         "Pr(>|t|)", "FDR")]
    results <- rbind(results, r_test)

    # Summarize
    sig_count <- sum(r_test$FDR < 0.05)
    cat("    Significant genes (FDR<0.05):", sig_count, "\n")
}

cat("\nStep 4: Summarizing results...\n")
# Overall summary
total_sig <- sum(results$FDR < 0.05)
high_fc_sig <- sum(results$FDR < 0.05 & abs(results$Estimate) > 1)

cat("  Total significant genes (FDR<0.05):", total_sig, "\n")
cat("  High FC genes (FDR<0.05, |log2FC|>1):", high_fc_sig, "\n")

# Top genes
cat("\nTop 10 genes by FDR:\n")
top_genes <- results[order(results$FDR), ]
print(head(top_genes[, c("Gene", "Subset", "Estimate", "FDR")], 10))

# Export results
cat("\nStep 5: Exporting results...\n")
write.csv(results, file = "geomx_de_results.csv", row.names = FALSE)
cat("  Full results: geomx_de_results.csv\n")

sig_genes <- subset(results, FDR < 0.05 & abs(Estimate) > 1)
write.csv(sig_genes, file = "geomx_de_significant.csv", row.names = FALSE)
cat("  Significant genes: geomx_de_significant.csv\n")

cat("\nAnalysis complete!\n")
