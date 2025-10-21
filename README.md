# Nectargenomics

A comparative genomics analysis pipeline associated with the preprint: [Convergent and lineage-specific genomic changes shape adaptations in sugar-consuming birds.](https://www.biorxiv.org/content/10.1101/2024.08.30.610474v1)

[![DOI](https://zenodo.org/badge/260459714.svg)](https://doi.org/10.5281/zenodo.17403724)

## Overview

This repository contains genomic analysis workflows for investigating the molecular basis of nectarivory adaptation across multiple bird species. The analyses include:

- **Assembly_and_annotation**: Genome assembly and annotation
- **CNEE Analysis (2021)**: Identification and analysis of conserved non-coding elements across bird genomes
- **Convergent Substitution Analysis (csubst)**: Detection of convergent molecular evolution
- **Positive Selection Analysis (absrel)**: Adaptive branch-site random effects likelihood tests
- **Transcriptome Analysis**: Gene expression data for nectar-feeding birds

## Repository Structure

```
nectargenomics/
├── Assembly_and_annotation  # Genome assembly and annotation (submodule)
├── absrel/                  # Positive selection analysis (submodule)
├── csubst/                  # Convergent substitution analysis
│   ├── 00_inputs/           # Input sequences and species lists
│   └── 01_run_csubst/       # csubst pipeline execution
├── CNEEanalysis_2021/       # Conserved non-coding element analysis
│   ├── 00_inputs/           # Input data and species lists
│   ├── 01a_phylofit/        # Phylogenetic model fitting
│   ├── 01b_cnee_prep/       # CNEE identification and preparation
│   ├── 02_cnee_liftover/    # Coordinate conversion between assemblies
│   └── 03_phyloacc/         # PhyloAcc analysis for accelerated evolution
├── Transcriptome_for_NectarGenomics/  # Transcriptome data (submodule)
└── LICENSE                  # GNU GPL v3 license
```
