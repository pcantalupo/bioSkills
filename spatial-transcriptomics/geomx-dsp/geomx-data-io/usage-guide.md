# GeoMx Data I/O - Usage Guide

## Overview

This skill covers loading and accessing NanoString GeoMx Digital Spatial Profiler data using GeomxTools in R/Bioconductor.

## Prerequisites

```r
BiocManager::install(c('GeomxTools', 'NanoStringNCTools'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Load my GeoMx DCC files"
- "Read GeoMx protein data"

## Example Prompts

### Loading Data
> "Load GeoMx RNA data from DCC files"

> "Read GeoMx protein data with analyte='protein'"

> "Aggregate my GeoMx probes to target level"

### Accessing Data
> "Show me the sample metadata from my GeoMx data"

> "Get the count matrix from my GeoMx object"

> "Extract only endogenous probes"

## What the Agent Will Do

1. Locate DCC files, PKC files, and annotation file
2. Load data into NanoStringGeoMxSet object
3. Verify data structure and dimensions
4. Optionally aggregate probes to targets
5. Provide access to counts, metadata, and QC metrics

## Tips

- **File Organization** - Keep DCCs in one folder, PKCs in another, annotations as Excel/CSV
- **RNA vs Protein** - RNA is default; specify `analyte = "protein"` for protein data
- **Probe vs Target** - RNA data loads as probes; use `aggregateCounts()` for gene-level
- **Protein Auto-Aggregation** - Protein data is automatically aggregated to target level
- **Annotation File** - Must have Sample_ID column matching DCC file names
