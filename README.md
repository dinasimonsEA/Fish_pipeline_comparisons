# Fish pipeline comparisons

## Overview

This repository contains a set of bioinformatics workflows and analysis notebooks developed to compare alternative sequence inference and taxonomic assignment approaches for freshwater fish environmental DNA (eDNA) metabarcoding datasets. The project evaluates how different sample inference and taxonomic classification methods influence biodiversity assessments across multiple fish datasets. Four datasets are included: Marchamley, Loch Insh (Ring Trial), Windermere and simulated mock dataset. The analysis compares outputs from three sequence inference pipelines and five taxonomic assignment methods, enabling assessment of concordance between bioinformatic approaches and known species records.

## Repository workflow

The project is organised into three main stages:

1. **Bioinformatics pipeline processing**
2. **Data formatting and integration**
3. **Comparison of bioinformatics methods**

The workflow begins with raw sequence data, generates ASV or OTU tables using multiple pipelines, performs taxonomic assignment, combines all outputs into standardised datasets, and finally compares the impact of software choices on ecological conclusions.

## Notebooks

### 1. Bioinformatics pipeline testing

This notebook implements and evaluates sequence inference (MetaBEAT, DADA2 and VSEARCH) and taxonomic assignment (MetaBEAT, BLAST + LCA, RDP, VSEARCH + LCA, MAPseq) tools for fish metabarcoding data. The workflow processes four independent datasets, from raw reads to taxonomic assignment. 

Sequence processing and inference steps include:

- Removal of ambiguous bases (Ns)
- Primer trimming using **cutadapt**
- Read filtering
- Error modelling
- Dereplication
- Sequence denoising
- Read merging
- Chimera removal

Taxonomic assignment is performed using a custom reference database combining modified **MetaFishLib** (fish references) **MIDORI2** (non-fish references).

The notebook generates:
- ASV count tables
- ASV sequence tables
- Taxonomic assignments from multiple methods

## 2. Format processed data

This notebook combines outputs from all sequence processing and taxonomic assignment pipelines into standardised datasets for downstream analysis. The aim is to convert multiple heterogeneous outputs into:

- A tidy long-format master dataset
- A phyloseq object for ecological analysis

It also retrieves supporting taxonomic metadata from WoRMS (World Register of Marine Species).

### Input Files

For each combination of:

- 4 datasets
- 3 denoising methods

the following files are imported.

### ASV Information

Before taxonomic assignment:

1. ASV count table
2. ASV sequence table

### Taxonomic Assignment Outputs

1. MetaBEAT
2. RDP
3. MAPseq
4. BLAST (LCA condensed)
5. VSEARCH (LCA condensed)

### Data Scale

The workflow processes:

- 4 datasets
- 3 denoising methods
- 7 output files per combination

Resulting in:

**84 files (4 × 3 × 7)**

These are merged into a single unified dataset suitable for comparative analyses.

### Special Handling

MetaBEAT outputs differ in structure from the other pipelines and must be cleaned before integration.

### Outputs

- Long-format master data frame
- phyloseq object
- Standardised taxonomy and metadata tables


## 3. Compare bioinformatics tools

### Purpose

This notebook evaluates how bioinformatic processing choices affect biodiversity results.

Using the standardised datasets generated in Notebook 2, the workflow compares both denoising and taxonomic assignment approaches.

### Denoising Methods Compared

- DADA2
- VSEARCH
- MetaBEAT

### Taxonomic Assignment Methods Compared

- RDP
- BLAST LCA
- MAPseq
- MetaBEAT
- VSEARCH LCA

### Comparison Metrics

The notebook compares:

- ASVs
- Unique sequences
- Species-level assignments
- Genus-level assignments
- Family-level assignments

### Datasets Evaluated

Analyses are conducted using datasets with known fish communities based on historical survey information:

- Ring Trial
- Marchamley
- Windermere

### Ecological Validation

Results are compared against visual survey records to assess:

- Taxonomic concordance
- Species recovery
- Method-specific biases
- Consistency across datasets

# Repository Structure

```text
Fish_pipeline_comparisons/
│
├── 01_Bioinformatics_pipeline_testing.ipynb
├── 02_Format_processed_data.ipynb
├── 03_Compare_bioinformatics_tools.ipynb
│
├── data/
├── outputs/
├── scripts/
└── README.md
```

# Key Comparisons

The project evaluates:

| Category | Methods |
|-----------|----------|
| Denoising | DADA2, VSEARCH, MetaBEAT |
| Taxonomy Assignment | RDP, BLAST LCA, MAPseq, MetaBEAT, VSEARCH LCA |
| Datasets | Marchamley, Ring Trial, Windermere, Simulated |
| Taxonomic Levels | Species, Genus, Family |

# Outputs

The repository produces:

- ASV abundance tables
- ASV sequence tables
- Taxonomic assignment tables
- Long-format integrated datasets
- phyloseq objects
- Comparative analyses and visualisations
- Concordance assessments against historical fish records

# Running the Analysis

The notebooks should be run sequentially:

1. **Bioinformatics Pipeline Testing** – process raw sequence data and generate denoised ASVs and taxonomic assignments.
2. **Format Processed Data** – combine outputs from all pipelines into standardised long-format datasets and phyloseq objects.
3. **Compare Bioinformatics Tools** – evaluate differences between denoising methods and taxonomic assignment approaches.

Where required, scripts can also be run independently outside the notebook workflow.