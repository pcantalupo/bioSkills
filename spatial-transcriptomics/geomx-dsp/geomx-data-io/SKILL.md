---
name: bio-spatial-transcriptomics-geomx-data-io
description: Load NanoString GeoMx Digital Spatial Profiler (DSP) data from DCC files using GeomxTools. Read RNA and protein expression data, PKC probe configuration files, and sample annotations. Use when loading GeoMx DSP NGS or protein data.
tool_type: r
primary_tool: GeomxTools
---

# GeoMx Data I/O

Load and work with NanoString GeoMx Digital Spatial Profiler data.

## Required Libraries

```r
library(GeomxTools)
library(NanoStringNCTools)
```

## Installation

```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')
BiocManager::install('GeomxTools')
BiocManager::install('NanoStringNCTools')
```

## Load GeoMx RNA Data

```r
# Specify file locations
datadir <- "path/to/geomx/data"
DCCFiles <- dir(file.path(datadir, "dccs"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
PKCFiles <- unzip(zipfile = dir(file.path(datadir, "pkcs"), pattern = ".zip$",
                                full.names = TRUE, recursive = TRUE))
SampleAnnotationFile <- file.path(datadir, "annotation.xlsx")

# Load RNA data (default analyte)
demoData <- readNanoStringGeoMxSet(
    dccFiles = DCCFiles,
    pkcFiles = PKCFiles,
    phenoDataFile = SampleAnnotationFile,
    phenoDataSheet = "Template",
    phenoDataDccColName = "Sample_ID",
    protocolDataColNames = c("aoi", "roi"),
    experimentDataColNames = c("panel")
)

print(demoData)
# NanoStringGeoMxSet (storageMode: lockedEnvironment)
# assayData: 18642 features, 88 samples
# feature: Probe
# analyte: RNA
```

## Load GeoMx Protein Data

```r
# Load protein data (must specify analyte)
proteinData <- readNanoStringGeoMxSet(
    dccFiles = DCCFiles,
    pkcFiles = PKCFiles,
    phenoDataFile = SampleAnnotationFile,
    phenoDataSheet = "Annotations",
    phenoDataDccColName = "Sample_ID",
    protocolDataColNames = c("Tissue", "Segment_Type", "ROI.Size"),
    configFile = NULL,
    analyte = "protein",  # Specify protein analyte
    phenoDataColPrefix = "",
    experimentDataColNames = NULL
)

print(proteinData)
# feature: Target (protein data is auto-aggregated)
# analyte: Protein
```

## Access Expression Data

```r
# Get count matrix
counts <- assayData(demoData)[["exprs"]]
dim(counts)  # probes x samples

# View first few counts
counts[1:5, 1:3]
```

## Access Sample Metadata (phenoData)

```r
# Get phenotype data
sample_info <- pData(demoData)
head(sample_info)

# Access specific columns
slide_names <- pData(demoData)$`slide name`
regions <- pData(demoData)$region
segments <- pData(demoData)$segment
```

## Access Protocol Data (Sequencing QC)

```r
# Get protocol/sequencing metrics
protocol_info <- pData(protocolData(demoData))
head(protocol_info)

# Sequencing metrics
seq_metrics <- protocol_info[, c("Raw", "Trimmed", "Stitched",
                                 "Aligned", "DeduplicatedReads")]
```

## Access Feature Metadata (Probes/Targets)

```r
# Get probe/target information
probe_info <- fData(demoData)
head(probe_info)

# Key columns
probe_info[1:5, c("RTS_ID", "TargetName", "Module", "CodeClass")]

# Check feature type (Probe or Target)
featureType(demoData)  # "Probe" or "Target"
```

## Check Panel Information

```r
# Get PKC annotation
pkc_files <- annotation(demoData)
print(pkc_files)

# Get module names
modules <- gsub(".pkc", "", pkc_files)
print(modules)
```

## Aggregate Probes to Targets

```r
# For RNA data, aggregate probe-level counts to target-level
# (Protein data is automatically aggregated on load)
target_demoData <- aggregateCounts(demoData)

# Check feature type changed from Probe to Target
featureType(target_demoData)  # "Target"

# Dimensions change (fewer features)
dim(demoData)         # 8707 probes x 88 samples
dim(target_demoData)  # 1821 targets x 88 samples

# Expression is now gene-level
head(rownames(target_demoData))  # Gene symbols
```

## Subset Data

```r
# Subset by sample
demoData_subset <- demoData[, pData(demoData)$region == "glomerulus"]

# Subset by feature
target_genes <- c("PDCD1", "CD274", "IFNG", "CD8A")
demoData_genes <- target_demoData[target_genes, ]

# Subset using bracket notation
demoData_filtered <- demoData[1:1000, 1:50]
```

## Get Endogenous Probes Only

```r
# Remove control probes, keep only endogenous
endo_data <- endogenousSubset(demoData)
dim(endo_data)

# Can also combine with sample subsetting
endo_glom <- endogenousSubset(
    demoData,
    select = pData(demoData)$region == "glomerulus"
)
```

## Get Negative Control Probes

```r
# Get negative control probes
neg_data <- negativeControlSubset(demoData)
dim(neg_data)

# View negative probe names
fData(neg_data)$TargetName
```

## Access Multiple Assay Matrices

```r
# After normalization, multiple expression matrices are stored
assayDataElementNames(demoData)
# [1] "exprs" "q_norm" "neg_norm"

# Access specific normalized data
q_norm_counts <- assayDataElement(demoData, elt = "q_norm")
neg_norm_counts <- assayDataElement(demoData, elt = "neg_norm")
```

## Study Design Summary

```r
# Summarize sample distribution
library(dplyr)
pData(demoData) %>%
    group_by(region, segment) %>%
    summarise(n_samples = n())

# Check sample counts
table(pData(demoData)$region)
table(pData(demoData)$segment)
```

## Export Data

```r
# Save as RDS
saveRDS(target_demoData, file = "geomx_data.rds")

# Load back
target_demoData <- readRDS("geomx_data.rds")

# Export counts to CSV
write.csv(assayData(target_demoData)[["exprs"]],
          file = "geomx_counts.csv")

# Export sample metadata
write.csv(pData(target_demoData),
          file = "geomx_sample_info.csv")
```

## Common Data Checks

```r
# Check dimensions
dim(demoData)

# Check if samples match between expression and metadata
all(colnames(exprs(demoData)) == rownames(pData(demoData)))

# Check for missing values
sum(is.na(exprs(demoData)))

# Summary statistics
summary(colSums(exprs(demoData)))  # Total counts per sample
summary(rowSums(exprs(demoData)))  # Total counts per gene
```

## Related Skills

- geomx-qc - Quality control for segments and probes
- geomx-normalization - Normalize GeoMx data
- spatial-preprocessing - General spatial data QC
