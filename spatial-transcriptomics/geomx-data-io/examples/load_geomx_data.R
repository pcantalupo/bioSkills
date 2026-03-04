#!/usr/bin/env Rscript
# Load GeoMx DSP data from DCC files

library(GeomxTools)
library(NanoStringNCTools)

# Specify data directory
datadir <- "path/to/geomx/data"

# Locate files
DCCFiles <- dir(file.path(datadir, "dccs"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
PKCFiles <- unzip(zipfile = dir(file.path(datadir, "pkcs"), pattern = ".zip$",
                                full.names = TRUE, recursive = TRUE))
SampleAnnotationFile <- file.path(datadir, "annotation.xlsx")

cat("Found", length(DCCFiles), "DCC files\n")
cat("Found", length(PKCFiles), "PKC files\n")

# Load RNA data
demoData <- readNanoStringGeoMxSet(
    dccFiles = DCCFiles,
    pkcFiles = PKCFiles,
    phenoDataFile = SampleAnnotationFile,
    phenoDataSheet = "Template",
    phenoDataDccColName = "Sample_ID",
    protocolDataColNames = c("aoi", "roi"),
    experimentDataColNames = c("panel")
)

cat("\nLoaded GeoMx data:\n")
print(demoData)

cat("\nFeature type:", featureType(demoData), "\n")
cat("Analyte:", demoData@annotation[1], "\n")

# Aggregate probes to targets for RNA data
cat("\nAggregating probes to targets...\n")
target_demoData <- aggregateCounts(demoData)

cat("After aggregation:\n")
cat("  Features:", nrow(target_demoData), "\n")
cat("  Samples:", ncol(target_demoData), "\n")
cat("  Feature type:", featureType(target_demoData), "\n")

# Sample distribution
cat("\nSample distribution:\n")
print(table(pData(target_demoData)$region))

# Save processed data
saveRDS(target_demoData, file = "geomx_targets.rds")
cat("\nSaved to geomx_targets.rds\n")
