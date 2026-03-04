# GeoMx Differential Expression - Usage Guide

## Overview

This skill covers differential expression analysis for GeoMx data using linear mixed-effect models (LMM) that account for tissue subsampling.

## Prerequisites

```r
BiocManager::install(c('GeomxTools', 'lme4', 'lmerTest'))
install.packages(c('ggplot2', 'ggrepel', 'pheatmap'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Compare glomerulus vs tubule in my GeoMx data"
- "Find differentially expressed genes between DKD and normal"

## Example Prompts

### Within-Slide Comparisons
> "Compare two regions within the same tissue slide"

> "Find genes enriched in glomerulus vs tubule"

> "Run DE with random slope model"

### Between-Slide Comparisons
> "Compare disease vs healthy across different slides"

> "Find DKD-specific genes in glomeruli"

> "Run DE with random intercept model"

### Visualization
> "Create a volcano plot of my DE results"

> "Plot expression of NPHS1 across regions"

> "Make a heatmap of significant genes"

## What the Agent Will Do

1. Prepare log-transformed expression data
2. Set up appropriate model formula (random slope vs intercept)
3. Run linear mixed-effect model
4. Calculate FDR-adjusted p-values
5. Generate volcano plots and heatmaps
6. Export results to CSV

## Tips

- **Model Choice**: Use random slope for within-slide, random intercept for between-slide
- **Log Transform**: Always work with log2-transformed normalized data
- **Thresholds**: FDR < 0.05 and |log2FC| > 1 is typical
- **Multiple Testing**: FDR correction accounts for testing thousands of genes
- **Visualization**: Volcano plots show both significance and effect size
