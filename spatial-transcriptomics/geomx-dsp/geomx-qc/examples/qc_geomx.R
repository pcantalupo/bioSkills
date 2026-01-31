#!/usr/bin/env Rscript
# Quality control for GeoMx DSP data

library(GeomxTools)
library(ggplot2)

# Load data (assuming already loaded)
# demoData <- readRDS("geomx_data.rds")

cat("Starting QC workflow...\n")
cat("Input:", nrow(demoData), "features x", ncol(demoData), "samples\n\n")

# 1. Shift zero counts
cat("Step 1: Shifting zero counts to 1...\n")
demoData <- shiftCountsOne(demoData, useDALogic = TRUE)

# 2. Set segment QC flags
cat("Step 2: Setting segment QC flags...\n")
# minSegmentReads=1000: Minimum acceptable raw read count for segment quality
# percentAligned=75: 75% alignment ensures accurate quantification (lower suggests mapping issues)
# percentSaturation=50: 50% saturation means sufficient read depth to capture transcript diversity
QC_params <- list(
    minSegmentReads = 1000,
    percentTrimmed = 80,
    percentStitched = 80,
    percentAligned = 75,
    percentSaturation = 50,
    minNegativeCount = 1,
    maxNTCCount = 9000,
    minNuclei = 20,
    minArea = 1000
)

demoData <- setSegmentQCFlags(demoData, qcCutoffs = QC_params)

# Check results
QCResults <- protocolData(demoData)[["QCFlags"]]
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
    ifelse(sum(x) == 0L, "PASS", "WARNING")
})

cat("  PASS:", sum(QCResults$QCStatus == "PASS"), "\n")
cat("  WARNING:", sum(QCResults$QCStatus == "WARNING"), "\n\n")

# 3. Remove failed segments
cat("Step 3: Removing low-quality segments...\n")
passedQC <- demoData[, QCResults$QCStatus == "PASS"]
cat("  Retained:", ncol(passedQC), "segments\n\n")

# 4. Probe QC
cat("Step 4: Setting probe QC flags...\n")
passedQC <- setBioProbeQCFlags(
    passedQC,
    qcCutoffs = list(minProbeRatio = 0.1, percentFailGrubbs = 20),
    removeLocalOutliers = TRUE
)

ProbeQCResults <- fData(passedQC)[["QCFlags"]]
cat("  Global outliers:", sum(ProbeQCResults$GlobalGrubbsOutlier), "\n")
cat("  Low ratio:", sum(ProbeQCResults$LowProbeRatio), "\n\n")

# 5. Remove outlier probes
cat("Step 5: Removing outlier probes...\n")
passedQC <- passedQC[
    !fData(passedQC)[["QCFlags"]][, "LowProbeRatio"] &
    !fData(passedQC)[["QCFlags"]][, "GlobalGrubbsOutlier"], ]
cat("  Retained:", nrow(passedQC), "probes\n\n")

# 6. Aggregate to targets
cat("Step 6: Aggregating to target level...\n")
target_data <- aggregateCounts(passedQC)
cat("  Targets:", nrow(target_data), "\n\n")

cat("QC complete!\n")
cat("Final:", nrow(target_data), "targets x", ncol(target_data), "segments\n")

# Save QC-passed data
saveRDS(target_data, file = "geomx_qc_passed.rds")
cat("\nSaved to geomx_qc_passed.rds\n")
