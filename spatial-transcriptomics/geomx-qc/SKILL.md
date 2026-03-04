---
name: bio-spatial-transcriptomics-geomx-qc
description: Quality control for GeoMx DSP segments and probes. Flag low-quality ROI/AOI segments based on sequencing metrics, saturation, negative controls, and tissue characteristics. Remove outlier probes using Grubb's test. Use when performing QC on GeoMx data.
tool_type: r
primary_tool: GeomxTools
---

# GeoMx Quality Control

Quality control for GeoMx DSP segments (ROI/AOI) and probes.

## Required Libraries

```r
library(GeomxTools)
library(ggplot2)
```

## Shift Zero Counts

```r
# Shift counts of 0 to 1 for downstream log transformations
# useDALogic=TRUE uses Detection Above Background logic
demoData <- shiftCountsOne(demoData, useDALogic = TRUE)
```

## Set Segment QC Flags

```r
# Define QC cutoffs for segments
QC_params <- list(
    minSegmentReads = 1000,      # Minimum raw reads
    percentTrimmed = 80,         # Minimum % trimmed (80%)
    percentStitched = 80,        # Minimum % stitched (80%)
    percentAligned = 75,         # Minimum % aligned (80%)
    percentSaturation = 50,      # Minimum sequencing saturation (50%)
    minNegativeCount = 1,        # Minimum negative control counts (10)
    maxNTCCount = 9000,          # Maximum NTC counts (1000)
    minNuclei = 20,              # Minimum nuclei count (100)
    minArea = 1000               # Minimum segment area (5000)
)

# Apply QC flags
demoData <- setSegmentQCFlags(demoData, qcCutoffs = QC_params)
```

## View QC Results

```r
# Get QC flags
QCResults <- protocolData(demoData)[["QCFlags"]]

# Summarize QC results
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(
    Pass = colSums(!QCResults[, flag_columns]),
    Warning = colSums(QCResults[, flag_columns])
)

print(QC_Summary)
#                LowReads LowTrimmed LowStitched LowAligned ...
# Pass               231        235         235        229
# Warning              4          0           0          6

# Overall QC status
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
    ifelse(sum(x) == 0L, "PASS", "WARNING")
})

table(QCResults$QCStatus)
# PASS    WARNING
# 229          6
```

## Visualize Segment QC Metrics

```r
# Function to plot QC distributions
QC_histogram <- function(assay_data, annotation, fill_by, thr = NULL) {
    plt <- ggplot(assay_data,
                  aes_string(x = paste0("unlist(`", annotation, "`)"),
                             fill = fill_by)) +
        geom_histogram(bins = 50) +
        geom_vline(xintercept = thr, lty = "dashed", color = "black") +
        theme_bw() +
        facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
        labs(x = annotation, y = "Segments, #")
    plt
}

# Plot various QC metrics
QC_histogram(sData(demoData), "Trimmed (%)", "segment", 80)
QC_histogram(sData(demoData), "Aligned (%)", "segment", 75)
QC_histogram(sData(demoData), "Saturated (%)", "segment", 50)
QC_histogram(sData(demoData), "nuclei", "segment", 20)
```

## Calculate Negative GeoMean

```r
# Calculate geometric mean of negative probes per module
negativeGeoMeans <-
    esBy(negativeControlSubset(demoData),
         GROUP = "Module",
         FUN = function(x) {
             assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs")
         })

# Add to protocol data
protocolData(demoData)[["NegGeoMean"]] <- negativeGeoMeans

# Plot negative geomean
library(scales)
negCols <- paste0("NegGeoMean_", modules)
pData(demoData)[, negCols] <- sData(demoData)[["NegGeoMean"]]

QC_histogram(pData(demoData), negCols[1], "segment", 2)
```

## Remove Low-Quality Segments

```r
# Get segments that passed QC
QCResultsIndex <- which(apply(protocolData(demoData)[["QCFlags"]],
                              1L, function(x) sum(x) == 0L))

# Subset to passing segments
QCPassed <- demoData[, QCResultsIndex]

cat("Before QC:", ncol(demoData), "segments\n")
cat("After QC:", ncol(QCPassed), "segments\n")
# Before QC: 235 segments
# After QC: 229 segments
```

## Set Probe QC Flags

```r
# Set QC flags for probes (outlier detection)
# minProbeRatio=0.1: Probe geomean / target geomean threshold
# percentFailGrubbs=20: Max % segments where probe is Grubb's outlier
demoData <- setBioProbeQCFlags(
    demoData,
    qcCutoffs = list(minProbeRatio = 0.1, percentFailGrubbs = 20),
    removeLocalOutliers = TRUE  # Remove segment-specific outliers
)

# Get probe QC results
ProbeQCResults <- fData(demoData)[["QCFlags"]]

# Summarize probe QC
qc_df <- data.frame(
    Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0 &
                !ProbeQCResults$GlobalGrubbsOutlier)
)

print(qc_df)
#   Passed Global Local
#    18619      1    22
```

## Remove Outlier Probes

```r
# Remove probes that failed global QC
ProbeQCPassed <- subset(
    demoData,
    fData(demoData)[["QCFlags"]][, "LowProbeRatio"] == FALSE &
    fData(demoData)[["QCFlags"]][, "GlobalGrubbsOutlier"] == FALSE
)

cat("Before probe QC:", nrow(demoData), "probes\n")
cat("After probe QC:", nrow(ProbeQCPassed), "probes\n")
# Before probe QC: 18642 probes
# After probe QC: 18641 probes

demoData <- ProbeQCPassed
```

## Calculate Limit of Quantification (LOQ)

```r
# LOQ = geomean(NegProbe) * geoSD(NegProbe)^n
# Typically n=2 (2 standard deviations above background)
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module
LOQ <- data.frame(row.names = colnames(target_demoData))
for(module in modules) {
    vars <- paste0(c("NegGeoMean_", "NegGeoSD_"), module)
    if(all(vars %in% colnames(pData(target_demoData)))) {
        LOQ[, module] <- pmax(
            minLOQ,
            pData(target_demoData)[, vars[1]] *
                pData(target_demoData)[, vars[2]] ^ cutoff
        )
    }
}

pData(target_demoData)$LOQ <- LOQ
```

## Gene Detection Rate

```r
# Calculate which genes are detected above LOQ in each segment
LOQ_Mat <- c()
for(module in modules) {
    ind <- fData(target_demoData)$Module == module
    Mat_i <- t(esApply(target_demoData[ind, ], MARGIN = 1,
                       FUN = function(x) {
                           x > LOQ[, module]
                       }))
    LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}

# Ensure ordering
LOQ_Mat <- LOQ_Mat[fData(target_demoData)$TargetName, ]

# Calculate detection rate per segment
pData(target_demoData)$GenesDetected <- colSums(LOQ_Mat, na.rm = TRUE)
pData(target_demoData)$GeneDetectionRate <-
    pData(target_demoData)$GenesDetected / nrow(target_demoData)

# Categorize detection rates
pData(target_demoData)$DetectionThreshold <-
    cut(pData(target_demoData)$GeneDetectionRate,
        breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
        labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

# Plot detection rate distribution
ggplot(pData(target_demoData), aes(x = DetectionThreshold)) +
    geom_bar(aes(fill = region)) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
    theme_bw() +
    labs(x = "Gene Detection Rate", y = "Segments, #")
```

## Filter Segments by Detection Rate

```r
# Remove segments with <10% genes detected
target_demoData <- target_demoData[,
    pData(target_demoData)$GeneDetectionRate >= 0.1]

cat("Segments after detection filtering:", ncol(target_demoData), "\n")
```

## Filter Genes by Detection Rate

```r
# Calculate detection rate per gene
fData(target_demoData)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_demoData)$DetectionRate <-
    fData(target_demoData)$DetectedSegments / ncol(target_demoData)

# Keep genes detected in ≥10% of segments
# Also keep negative control probes
negativeProbefData <- subset(fData(target_demoData), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)

target_demoData <- target_demoData[
    fData(target_demoData)$DetectionRate >= 0.1 |
    fData(target_demoData)$TargetName %in% neg_probes, ]

cat("Genes after detection filtering:", nrow(target_demoData), "\n")
```

## Complete QC Workflow

```r
library(GeomxTools)

# Load data
demoData <- readNanoStringGeoMxSet(...)

# 1. Shift zero counts
demoData <- shiftCountsOne(demoData, useDALogic = TRUE)

# 2. Segment QC
demoData <- setSegmentQCFlags(demoData, qcCutoffs = QC_params)
QCResults <- protocolData(demoData)[["QCFlags"]]
passedQC <- demoData[, rowSums(QCResults) == 0]

# 3. Probe QC
passedQC <- setBioProbeQCFlags(passedQC)
passedQC <- passedQC[
    !fData(passedQC)[["QCFlags"]][, "LowProbeRatio"] &
    !fData(passedQC)[["QCFlags"]][, "GlobalGrubbsOutlier"], ]

# 4. Aggregate to targets
target_data <- aggregateCounts(passedQC)

# 5. Calculate LOQ and detection rates
# (see above sections)

# 6. Filter segments and genes by detection
target_data <- target_data[, pData(target_data)$GeneDetectionRate >= 0.1]
target_data <- target_data[fData(target_data)$DetectionRate >= 0.1, ]

cat("Final dataset:", nrow(target_data), "genes x",
    ncol(target_data), "segments\n")
```

## Related Skills

- geomx-data-io - Load GeoMx data
- geomx-normalization - Normalize after QC
- spatial-preprocessing - General spatial QC
