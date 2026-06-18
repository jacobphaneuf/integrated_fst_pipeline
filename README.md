# Integrated FST Pipeline

## Overview

This repository contains an integrated R pipeline for fecal source tracking (FST) in surface waters using sequencing, spatial, and digital PCR (dPCR) data.

Upstream processing was performed externally using Kraken2, Bracken, metaSourceTracker2 (mST2), and SourceApp. These tools are not implemented within this repository but generate inputs used by the pipeline.

## Features

* Sequencing normalization using:

  * CSS (Cumulative Sum Scaling via metagenomeSeq)
  * TMM (Trimmed Mean by M-Values via edgeR)
  * MED (Median Ratio Normalization via DESeq2)
  * VST (Variance Stabilizing Transformation via DESeq2)

* Beta diversity analyses:

  * Bray-Curtis dissimilarity
  * ANOSIM significance testing
  * Pairwise ANOSIM
  * NMDS visualization

* Source tracking tools:

  * Running FEAST in triplicate 
  * Preparing data for MetaSourceTracker2 (mST2)
  * Comparing tool outputs to SourceApp results

* Correlation analyses (Spearman's rank correlations with false discovery rate correction)

* Differential abundance testing (ANCOM-BC2)

* Publication-quality figure generation

## Required Input Files

The pipeline requires the following input data:

* Bracken Outputs (".br" Files)

* Metadata (User-Generated)

* Additional Data Files (Marker/Spatial Datasets)

## Pipeline Workflow

### Step 1

Normalization testing, beta diversity analyses, and NMDS generation.

### Step 2

Running FEAST and generation of mST2 input tables.

### Step 3

Correlation analyses between source tracking outputs and dPCR markers.

### Step 4

ANCOM-BC2 differential abundance analyses.

### Step 5

Assembly of publication-ready NMDS plots and heatmaps.

### Step 6

Visualization of SourceApp MAGs by state.

## Outputs

The pipeline generates:

* Normalized OTU tables
* ANOSIM statistics
* NMDS plots
* FEAST outputs
* Correlation statistics
* ANCOM-BC2 results
* Publication-quality TIFF figures

## Software Requirements

* R (version 4.0 or later)
* Required R packages listed within the script

## Citation

If you use this pipeline, please cite:

Phaneuf J.R. et al. (manuscript submitted, under review).
A citable DOI will be added upon publication.

## Software References

This pipeline implements or utilizes previously described methods and software packages from the following publications:

* **edgeR** — Chen Y, Chen L, Lun ATL, Baldoni PL, Smyth GK. 2025. *edgeR v4: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets*. *Nucleic Acids Research* 53(2):gkaf018. https://doi.org/10.1093/nar/gkaf018

* **ANCOM-BC2** — Lin H, Peddada SD. 2024. *Multigroup analysis of compositions of microbiomes with covariate adjustments and repeated measures*. *Nature Methods* 21(1):83–91. https://doi.org/10.1038/s41592-023-02092-7

* **DESeq2** — Love MI, Huber W, Anders S. 2014. *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2*. *Genome Biology* 15(12):550. https://doi.org/10.1186/s13059-014-0550-8

* **vegan** — Oksanen J, Simpson GL, Blanchet FG, et al. 2026. *vegan: Community Ecology Package*. R package version 2.7-3. https://doi.org/10.32614/CRAN.package.vegan

* **metagenomeSeq (CSS normalization)** — Paulson JN, Stine OC, Bravo HC, Pop M. 2013. *Differential abundance analysis for microbial marker-gene surveys*. *Nature Methods* 10:1200–1202. https://doi.org/10.1038/nmeth.2658

* **FEAST** — Shenhav L, Thompson M, Joseph TA, et al. 2019. *FEAST: fast expectation-maximization for microbial source tracking*. *Nature Methods* 16:627–632. https://doi.org/10.1038/s41592-019-0431-x
