# geomx-dsp

## Overview

Analyze NanoString GeoMx Digital Spatial Profiler (DSP) data for spatially-resolved transcriptomics and proteomics using GeomxTools in R/Bioconductor.

**Tool type:** r | **Primary tools:** GeomxTools, NanoStringNCTools

## Skills

| Skill | Description |
|-------|-------------|
| geomx-data-io | Load DCC files, PKC configuration, and sample annotations |
| geomx-qc | Quality control for segments (ROI/AOI) and probes |
| geomx-normalization | Q3, housekeeping, or background normalization |
| geomx-differential-expression | Linear mixed-effect models for DE analysis |

## Example Prompts

- "Load my GeoMx DCC files with protein data"
- "Run QC on GeoMx segments and remove low-quality ROIs"
- "Normalize my GeoMx RNA data using Q3 method"
- "Compare glomerulus vs tubule expression in GeoMx data"
- "Find genes differentially expressed between DKD and healthy kidneys"
- "Create a volcano plot of my GeoMx DE results"
- "Filter segments with less than 10% gene detection"
- "Calculate LOQ for my GeoMx data"

## Requirements

```r
BiocManager::install(c('GeomxTools', 'NanoStringNCTools', 'lme4', 'lmerTest'))
install.packages(c('ggplot2', 'ggrepel', 'pheatmap', 'reshape2', 'dplyr'))
```

## Workflow Overview

```
1. Data Loading (geomx-data-io)
   └─> Load DCC files, PKC, annotations
   └─> Aggregate probes to targets (RNA only)

2. Quality Control (geomx-qc)
   └─> Segment QC (sequencing metrics, saturation)
   └─> Probe QC (outlier detection)
   └─> LOQ calculation
   └─> Gene detection filtering

3. Normalization (geomx-normalization)
   └─> Q3 for RNA (recommended)
   └─> Housekeeping or background for protein
   └─> Log transformation

4. Differential Expression (geomx-differential-expression)
   └─> Linear mixed-effect models
   └─> Account for tissue subsampling
   └─> Volcano plots, heatmaps
```

## GeoMx vs Other Spatial Platforms

| Feature | GeoMx DSP | Visium | Xenium |
|---------|-----------|--------|--------|
| Resolution | ROI-based (user-defined) | Spot-based (55µm) | Single-cell |
| Throughput | NGS readout | NGS readout | In-situ imaging |
| Content | WTA or panels | WTA | Panels |
| Segments | Multiple per slide | Fixed grid | Subcellular |
| Analysis | GeomxTools (R) | Scanpy/Seurat | Squidpy/SpatialData |

## Key Differences from Standard RNA-seq

- **ROI/AOI Structure**: Multiple regions per tissue → Use mixed models
- **Probe-Level Data**: RNA data loads as probes → Aggregate to targets
- **Detection-Based Filtering**: Use LOQ rather than count thresholds
- **Normalization**: Q3 preferred over TMM/DESeq2 size factors
- **Protein Data**: Auto-aggregated, use HK or background normalization

## Related Skills

- **spatial-transcriptomics** - Visium, Xenium, other platforms
- **single-cell** - Non-spatial scRNA-seq analysis
- **differential-expression** - Alternative DE methods (DESeq2, edgeR)
- **data-visualization** - Volcano plots, heatmaps, QC plots
