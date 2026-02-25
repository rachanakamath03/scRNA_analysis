## Single-Cell RNA-seq Analysis of PBMCs in Systemic Lupus Erythematosus (SLE)

## Project Overview
This repository contains a complete single-cell RNA-seq (scRNA-seq) analysis pipeline performed using Seurat (R) to explore immune cell heterogeneity in Peripheral Blood Mononuclear Cells (PBMCs) from:
2 SLE patients
1 Healthy control
The objective was to perform an end-to-end analysis including quality control, clustering, marker identification, and differential expression to characterize immune alterations in Systemic Lupus Erythematosus.

## Dataset Information
GEO Accession: GSE162577
Platform: 10x Genomics Chromium
Technology: Single-cell RNA sequencing
Samples:
SLE1
SLE2
CTRL
The dataset investigates immune cell heterogeneity in SLE PBMCs.

## Analysis Workflow
The following steps were performed:
1. Quality Control
Removal of low-quality cells
Filtering thresholds:
nFeature_RNA > 200
nFeature_RNA < 6000
percent.mt < 15%
2. Normalization
LogNormalize (scale factor = 10,000)
3. Feature Selection
2,000 highly variable genes (vst method)
4. Dimensionality Reduction
PCA
Elbow plot diagnostics
UMAP embedding
5. Clustering
Graph-based clustering
Resolution = 0.5
PCs used: 1–12
6. Cell Type Annotation
Marker-based annotation using canonical immune markers:
CD3D, CD3E → T cells
S100A8 → Inflammatory monocytes
NCR1 → NK cells
MS4A1 → B cells
7. Differential Expression
Comparison:
SLE vs Control
Parameters:
logfc.threshold = 0.25
min.pct = 0.1
Key Findings
Identification of 18 transcriptionally distinct immune clusters.
Expanded inflammatory monocyte signatures in SLE (e.g., S100A8 upregulation).
Altered T cell gene expression (CD3D/CD3E downregulation).
Presence of rare clusters potentially representing interferon-responsive subsets.

## Requirements
R version ≥ 4.2
Required Packages
Seurat
dplyr
ggplot2
