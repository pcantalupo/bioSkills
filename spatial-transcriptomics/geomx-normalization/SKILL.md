---
name: bio-spatial-transcriptomics-geomx-normalization
description: Normalize GeoMx DSP expression data using quantile (Q3), housekeeping gene, or background (negative probe) normalization methods. For RNA use Q3; for protein use housekeeping or background. Use when normalizing GeoMx expression data.
tool_type: r
primary_tool: GeomxTools
---

# GeoMx Normalization

Normalize GeoMx DSP expression data for downstream analysis.

## Required Libraries

```r
library(GeomxTools)
library(ggplot2)
library(reshape2)
```

## Normalization Methods Overview

| Method | Best For | Description |
|--------|----------|-------------|
| **Q3 (Quantile)** | RNA-seq | 75th percentile normalization (recommended for RNA) |
| **Background** | Protein, RNA | Normalize by negative control geometric mean |
| **Housekeeping** | Protein | Normalize by housekeeping gene geometric mean |

## Q3 Normalization (RNA - Recommended)

```r
# Quantile normalization to 75th percentile
# desiredQuantile=0.75: Upper quartile captures biological signal while being robust to outliers
target_demoData <- normalize(
    target_demoData,
    norm_method = "quant",
    desiredQuantile = 0.75,
    toElt = "q_norm"
)

# Access normalized data
q_norm_counts <- assayDataElement(target_demoData, elt = "q_norm")
```

## Background Normalization

```r
# Normalize by negative control probes
# Appropriate for both RNA and protein
target_demoData <- normalize(
    target_demoData,
    norm_method = "neg",
    fromElt = "exprs",  # Source data (raw counts)
    toElt = "neg_norm"  # Destination slot
)

# Access background-normalized data
neg_norm_counts <- assayDataElement(target_demoData, elt = "neg_norm")
```

## Housekeeping Gene Normalization

```r
# Normalize by housekeeping gene geometric mean
# Recommended for protein data
target_demoData <- normalize(
    target_demoData,
    norm_method = "hk",
    fromElt = "exprs",
    toElt = "hk_norm"
)

# Access housekeeping-normalized data
hk_norm_counts <- assayDataElement(target_demoData, elt = "hk_norm")
```

## Check Available Normalization Matrices

```r
# View all expression matrices in object
assayDataElementNames(target_demoData)
# [1] "exprs"    "q_norm"   "neg_norm" "hk_norm"

# Raw counts are always in "exprs"
# Normalized counts in respective slots
```

## Visualize Pre-Normalization QC

```r
# Compare Q3 values to negative probe geometric mean
# Good separation indicates stable Q3 normalization
library(reshape2)

ann_of_interest <- "region"
neg_probes <- subset(fData(target_demoData), CodeClass == "Negative")$TargetName

Stat_data <- data.frame(
    row.names = colnames(exprs(target_demoData)),
    Segment = colnames(exprs(target_demoData)),
    Annotation = pData(target_demoData)[, ann_of_interest],
    Q3 = unlist(apply(exprs(target_demoData), 2, quantile, 0.75, na.rm = TRUE)),
    NegProbe = exprs(target_demoData)[neg_probes, ]
)

Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"))

# Plot distributions
ggplot(Stat_data_m, aes(x = Value, fill = Statistic)) +
    geom_histogram(bins = 40) +
    theme_bw() +
    scale_x_continuous(trans = "log2") +
    facet_wrap(~Annotation, nrow = 1) +
    labs(x = "Counts", y = "Segments, #")

# Q3 vs NegProbe scatter
ggplot(Stat_data, aes(x = NegProbe, y = Q3, color = Annotation)) +
    geom_abline(intercept = 0, slope = 1, lty = "dashed") +
    geom_point() +
    theme_bw() +
    scale_x_continuous(trans = "log2") +
    scale_y_continuous(trans = "log2") +
    labs(x = "Negative Probe GeoMean", y = "Q3 Value")
```

## Visualize Post-Normalization

```r
# Boxplots of raw vs normalized data
# First 10 segments
par(mfrow = c(1, 3))

# Raw counts
boxplot(exprs(target_demoData)[, 1:10],
        col = "#9EDAE5",
        main = "Raw Counts",
        log = "y",
        names = 1:10,
        ylab = "Counts, Raw")

# Q3 normalized
boxplot(assayDataElement(target_demoData[, 1:10], elt = "q_norm"),
        col = "#2CA02C",
        main = "Q3 Norm Counts",
        log = "y",
        names = 1:10,
        ylab = "Counts, Q3 Normalized")

# Background normalized
boxplot(assayDataElement(target_demoData[, 1:10], elt = "neg_norm"),
        col = "#FF7F0E",
        main = "Neg Norm Counts",
        log = "y",
        names = 1:10,
        ylab = "Counts, Neg. Normalized")
```

## Protein-Specific: Check Housekeeping Genes

```r
# For protein data, check available HK genes
hk_names <- hkNames(proteinData)
print(hk_names)
# [1] "Histone H3" "GAPDH"      "S6"

# Get IgG negative control names
igg_names <- iggNames(proteinData)
print(igg_names)
# [1] "Ms IgG1"  "Ms IgG2a" "Rb IgG"
```

## Protein-Specific: QC Signal Quality

```r
# For protein: check signal vs background
# Visualize protein targets relative to IgG controls
fig <- qcProteinSignal(object = proteinData, neg.names = igg_names)
fig()

# Get protein order by signal strength
proteinOrder <- qcProteinSignalNames(object = proteinData, neg.names = igg_names)
head(proteinOrder, 20)
```

## Protein-Specific: Normalization Concordance

```r
# Calculate different normalization factors
normfactors <- computeNormalizationFactors(
    object = proteinData,
    area = "AOI.Size.um2",
    nuclei = "Nuclei.Counts"
)

# Plot concordance between normalization methods
plotNormFactorConcordance(
    object = proteinData,
    plotFactor = "Tissue",
    normfactors = normfactors
)
```

## Multiple Normalization Strategy

```r
# For comprehensive analysis, generate all three
target_demoData <- normalize(target_demoData, norm_method = "quant",
                             desiredQuantile = 0.75, toElt = "q_norm")
target_demoData <- normalize(target_demoData, norm_method = "neg",
                             fromElt = "exprs", toElt = "neg_norm")
target_demoData <- normalize(target_demoData, norm_method = "hk",
                             fromElt = "exprs", toElt = "hk_norm")

# Compare results
assayDataElementNames(target_demoData)
# [1] "exprs"    "q_norm"   "neg_norm" "hk_norm"

# Use Q3 for RNA downstream analyses
# Use HK or neg for protein downstream analyses
```

## Log-Transform for Visualization

```r
# Create log2-transformed version for visualization
# Add 1 to avoid log(0)
assayDataElement(object = target_demoData, elt = "log_q") <-
    assayDataApply(target_demoData, 2, FUN = log, base = 2, elt = "q_norm")

# Now log_q contains log2(Q3-normalized counts)
log_counts <- assayDataElement(target_demoData, elt = "log_q")
```

## Normalization Recommendations

### RNA Data
- **Primary**: Q3 (quantile) normalization
- **Alternative**: Background if negative probe counts are stable
- **Not recommended**: Housekeeping (unless validated HK genes available)

### Protein Data
- **Primary**: Housekeeping gene normalization
- **Alternative**: Background (IgG) normalization
- **Not recommended**: Q3 for protein

## Complete Normalization Workflow

```r
library(GeomxTools)

# Assuming QC-passed data
target_demoData <- readRDS("geomx_qc_passed.rds")

# RNA normalization
if (annotation(target_demoData)[1] == "RNA") {
    # Q3 normalization for RNA
    target_demoData <- normalize(
        target_demoData,
        norm_method = "quant",
        desiredQuantile = 0.75,
        toElt = "q_norm"
    )

    # Log transform for visualization
    assayDataElement(target_demoData, elt = "log_q") <-
        assayDataApply(target_demoData, 2, FUN = log, base = 2, elt = "q_norm")
}

# Protein normalization
if (annotation(target_demoData)[1] == "Protein") {
    # HK normalization for protein
    target_demoData <- normalize(
        target_demoData,
        norm_method = "hk",
        fromElt = "exprs",
        toElt = "hk_norm"
    )

    # Alternative: background normalization
    target_demoData <- normalize(
        target_demoData,
        norm_method = "neg",
        fromElt = "exprs",
        toElt = "neg_norm"
    )
}

# Save normalized data
saveRDS(target_demoData, file = "geomx_normalized.rds")
```

## Related Skills

- geomx-qc - QC before normalization
- geomx-differential-expression - DE analysis after normalization
- differential-expression/deseq2-basics - Alternative DE methods
