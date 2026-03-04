#!/usr/bin/env Rscript
# Normalize GeoMx DSP data

library(GeomxTools)
library(ggplot2)

# Load QC-passed data
# target_demoData <- readRDS("geomx_qc_passed.rds")

cat("Normalization workflow\n")
cat("Data:", nrow(target_demoData), "genes x", ncol(target_demoData), "segments\n\n")

# Check analyte type
analyte_type <- target_demoData@annotation[1]
cat("Analyte type:", analyte_type, "\n\n")

if (grepl("RNA", analyte_type, ignore.case = TRUE)) {
    cat("Applying Q3 normalization for RNA...\n")

    # Q3 normalization
    # desiredQuantile=0.75: 75th percentile is robust to outliers while capturing biological signal
    target_demoData <- normalize(
        target_demoData,
        norm_method = "quant",
        desiredQuantile = 0.75,
        toElt = "q_norm"
    )

    cat("  Q3 normalization complete\n")

    # Create log-transformed version
    assayDataElement(target_demoData, elt = "log_q") <-
        assayDataApply(target_demoData, 2, FUN = log, base = 2, elt = "q_norm")

    cat("  Log2 transformation complete\n\n")

    # Show summary statistics
    q3_values <- apply(assayDataElement(target_demoData, elt = "q_norm"),
                      2, quantile, 0.75, na.rm = TRUE)
    cat("  Q3 summary:\n")
    print(summary(q3_values))

} else {
    cat("Applying housekeeping normalization for protein...\n")

    # Housekeeping normalization for protein
    target_demoData <- normalize(
        target_demoData,
        norm_method = "hk",
        fromElt = "exprs",
        toElt = "hk_norm"
    )

    cat("  Housekeeping normalization complete\n")

    # Also apply background normalization
    cat("Applying background normalization...\n")
    target_demoData <- normalize(
        target_demoData,
        norm_method = "neg",
        fromElt = "exprs",
        toElt = "neg_norm"
    )

    cat("  Background normalization complete\n\n")
}

# Show available expression matrices
cat("Available expression matrices:\n")
print(assayDataElementNames(target_demoData))

cat("\nNormalization complete!\n")

# Save normalized data
saveRDS(target_demoData, file = "geomx_normalized.rds")
cat("\nSaved to geomx_normalized.rds\n")
