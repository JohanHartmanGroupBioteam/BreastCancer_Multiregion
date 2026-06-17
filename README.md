# BreastCancer_Multiregion

## Overview
This repository contains the main analysis pipelines for bulk DNA and spatial transcriptomics data, combining Visium spatial transcriptomics, Xenium in situ sequencing, and bulk DNA CNV analysis.


## Analysis Modules

### 1. Visium-Niche Analysis
Spatial analysis of tissue niches using 10x Visium data.
- **Key features:**
  - Spatial clustering and niche identification
  - Ligand-receptor interaction analysis
  - Spatial trajectory inference
- **Dependencies:** Seurat

### 2. Xenium-CNV Analysis
Copy number variation analysis from Xenium in situ sequencing data.
- **Key features:**
  - Single-cell CNV detection
  - Spatial CNV heterogeneity mapping
  - Clonal evolution tracking
- **Dependencies:** pandas, scanpy, infercnvpy

### 3. Bulk DNA-CNV (Read Counting Workflow)
Traditional read-depth based CNV detection from bulk sequencing.
- **Key features:**
  - Read counting and normalization
  - Segmentation algorithms
  - Visualization of CNV profiles
- **Dependencies:** FACETS, WhatshAP, SHAPEIT



