# spatial-transcriptomics

## Overview

Analyze spatial transcriptomics data from Visium, Xenium, MERFISH, GeoMx DSP, and other platforms using Squidpy, SpatialData, and GeomxTools.

**Tool type:** python, r | **Primary tools:** Squidpy, SpatialData, Scanpy, scimap, GeomxTools

## Skills

| Skill | Description |
|-------|-------------|
| spatial-data-io | Load spatial data from Visium, Xenium, Slide-seq, MERFISH |
| spatial-preprocessing | QC, normalization, and feature selection for spatial data |
| spatial-neighbors | Build spatial neighbor graphs and compute connectivity |
| spatial-statistics | Moran's I, spatial autocorrelation, co-occurrence, enrichment |
| spatial-domains | Identify spatial domains and tissue regions |
| image-analysis | Process and analyze tissue images with Squidpy |
| spatial-visualization | Static and interactive visualization of spatial data |
| spatial-communication | Ligand-receptor analysis and cell-cell interactions |
| spatial-deconvolution | Estimate cell type composition per spot |
| spatial-multiomics | Analyze high-resolution platforms (Slide-seq, Stereo-seq, Visium HD) |
| spatial-proteomics | Analyze CODEX, IMC, MIBI spatial proteomics data |
| geomx-data-io | Load GeoMx DCC files, PKC configuration, and sample annotations |
| geomx-qc | Quality control for GeoMx segments (ROI/AOI) and probes |
| geomx-normalization | Q3, housekeeping, or background normalization for GeoMx |
| geomx-differential-expression | Linear mixed-effect models for GeoMx DE analysis |

## Example Prompts

- "Load my Visium data"
- "Read this Xenium output folder"
- "Run QC on my spatial data"
- "Normalize my spatial transcriptomics data"
- "Build a spatial neighbor graph with 6 neighbors"
- "Calculate Moran's I for this gene"
- "Find spatially variable genes"
- "Run co-occurrence analysis"
- "Identify spatial domains in my tissue"
- "Segment cells from the H&E image"
- "Plot gene expression on the tissue"
- "Show clusters overlaid on the image"
- "Run ligand-receptor analysis"
- "Deconvolve my Visium data with cell2location"
- "Analyze my Slide-seq data"
- "Process Stereo-seq at bin level"
- "Work with Visium HD subcellular resolution"
- "Analyze my CODEX spatial proteomics data"
- "Find spatial interactions between cell types in IMC data"
- "Load my GeoMx DCC files"
- "Run QC on GeoMx segments and filter low saturation ROIs"
- "Normalize GeoMx RNA data with Q3 method"
- "Compare glomerulus vs tubule in my GeoMx kidney data"
- "Find differentially expressed genes between disease and healthy in GeoMx"

## Requirements

### Python (Visium, Xenium, MERFISH, etc.)
```bash
pip install squidpy spatialdata spatialdata-io scanpy anndata scimap
```

### R/Bioconductor (GeoMx DSP)
```r
BiocManager::install(c('GeomxTools', 'NanoStringNCTools'))
```

## Related Skills

- **single-cell** - Non-spatial scRNA-seq analysis
- **differential-expression** - DE between spatial regions
- **data-visualization** - Visualization of spatial patterns
