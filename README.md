# Mannix DGE & Single-Cell Pipeline

This repository contains the bulk RNA-seq differential gene expression (DGE) and single-cell/single-nucleus RNA-seq analysis pipelines used for the Mannix project.

## 1. Differential Gene Expression

The DGE workflow performs bulk RNA-seq differential expression analysis.

Main analyses include:

- DESeq2 differential gene expression
- PCA
- volcano plots
- gene-level visualization
- pathway and enrichment analyses

Run the DGE pipeline with:

```bash
snakemake --use-conda --profile ../profile_slurm/
```

---

## 2. Single-cell / Single-nucleus RNA-seq

The single-cell workflow analyzes the GSE304446 dataset using Seurat.

Main steps include:

- download GEO data
- create the Seurat object
- normalization
- variable feature selection
- scaling
- PCA
- UMAP
- gene-expression visualization
- pseudobulk aggregation by sample and cell type
- differential expression with DESeq2

The Seurat object is generated from:

- `GSE304446_matrix.mtx.gz`
- `GSE304446_features.tsv.gz`
- `GSE304446_barcodes.tsv.gz`
- `GSE304446_cell_metadata.tsv.gz`

### Download the single-cell data

The GEO files are downloaded locally:

```bash
snakemake download_single_cell --profile ../profile_local/
```

### Run the single-cell analysis

```bash
snakemake --use-conda --profile ../profile_slurm/
```

The workflow follows approximately:

```text
GEO data
   |
   v
Create Seurat object
   |
   +------> UMAP
   |
   +------> Gene-expression / pseudobulk analysis
                 |
                 v
               DESeq2
                 |
                 v
                :D
```

