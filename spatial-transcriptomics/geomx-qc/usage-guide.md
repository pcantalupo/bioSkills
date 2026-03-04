# GeoMx QC - Usage Guide

## Overview

This skill covers quality control for GeoMx DSP segments and probes, including sequencing QC, probe outlier detection, and LOQ calculation.

## Prerequisites

```r
BiocManager::install(c('GeomxTools', 'ggplot2'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Run QC on my GeoMx segments"
- "Filter low-quality ROIs from my GeoMx data"

## Example Prompts

### Segment QC
> "Flag segments with low sequencing quality"

> "Remove segments with less than 50% saturation"

> "Calculate LOQ for each segment"

### Probe QC
> "Find outlier probes using Grubb's test"

> "Remove probes with low probe ratio"

### Detection
> "Calculate gene detection rate per segment"

> "Filter genes detected in less than 10% of segments"

## What the Agent Will Do

1. Shift zero counts to 1
2. Set segment QC flags based on sequencing metrics
3. Visualize QC distributions
4. Remove low-quality segments
5. Detect probe outliers using Grubb's test
6. Calculate limit of quantification (LOQ)
7. Filter segments and genes by detection rate

## Tips

- **Study-Specific Thresholds** - Adjust QC cutoffs based on tissue type and experimental design
- **Saturation** - GeoMx typically requires ≥50% saturation; <50% may need re-sequencing
- **Nuclei Count** - Threshold depends on segment size; smaller AOIs may have fewer nuclei
- **Detection Rate** - 5-10% is typical for segment filtering; 10-20% for gene filtering
- **LOQ Calculation** - 2 geometric SDs above negative control background is standard
