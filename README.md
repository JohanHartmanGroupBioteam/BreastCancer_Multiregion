# BreastCancer_Multiregion

## Overview
This repository contains the main analysis pipelines for bulk DNA and spatial transcriptomics data, combining Visium spatial transcriptomics, Xenium in situ sequencing, and bulk DNA CNV analysis.


## Analysis Modules

### 1. Visium-Niche Analysis
Spatial analysis of tissue niches using 10x Visium data.
- **Key features:**
  - Preprocessing and deconvolution
  - Spatial clustering and niche identification
  - Spot-level intrinsic subtyping
- **Dependencies:** Seurat, BreastSubtypeR, nicheDE

### 2. Xenium-CNV Analysis
Copy number variation analysis from Xenium in situ sequencing data.
- **Key features:**
  - Smooth and normalization
  - Spatial CNV heterogeneity mapping
- **Dependencies:** pandas, scanpy, infercnvpy

### 3. Bulk DNA-CNV (Read Counting Workflow)
Read-depth based CNV detection from bulk sequencing, with haplotype-specific read counting for improved sensitivity and resolution.
- **Key features:**
  - haplotype-specific phasing using germline variants
  - Allele-specific read counting and binning
- **Dependencies:** FACETS, WhatshAP, SHAPEIT



