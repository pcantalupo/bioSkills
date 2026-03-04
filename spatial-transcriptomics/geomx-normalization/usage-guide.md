# GeoMx Normalization - Usage Guide

## Overview

This skill covers normalization of GeoMx DSP expression data using Q3, housekeeping, or background methods.

## Prerequisites

```r
BiocManager::install(c('GeomxTools', 'ggplot2'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Normalize my GeoMx RNA data using Q3"
- "Apply housekeeping normalization to protein data"

## Example Prompts

### RNA Normalization
> "Normalize my GeoMx data with Q3 method"

> "Use quantile normalization at 75th percentile"

> "Check if Q3 and negative probes are well separated"

### Protein Normalization
> "Normalize protein data with housekeeping genes"

> "Apply IgG background normalization"

> "Check protein signal quality before normalization"

## What the Agent Will Do

1. Check analyte type (RNA or protein)
2. Visualize Q3 vs negative probe separation
3. Apply appropriate normalization method
4. Create log-transformed data for visualization
5. Generate normalization QC plots

## Tips

- **RNA**: Use Q3 normalization (recommended)
- **Protein**: Use housekeeping or background normalization
- **QC First**: Always run normalization after segment/probe QC
- **Separation**: Ensure Q3 >> negative probes before Q3 normalization
- **Multiple Methods**: Can apply multiple normalizations and compare
