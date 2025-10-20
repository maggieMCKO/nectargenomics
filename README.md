# Nectargenomics

A comprehensive genomics analysis pipeline for studying the evolution of nectarivory in birds, with a focus on conserved non-coding elements (CNEEs), convergent substitutions, and positive selection.

## Overview

This repository contains genomic analysis workflows for investigating the molecular basis of nectarivory adaptation across multiple bird species. The analyses include:

- **CNEE Analysis (2021)**: Identification and analysis of conserved non-coding elements across bird genomes
- **Convergent Substitution Analysis (csubst)**: Detection of convergent molecular evolution
- **Positive Selection Analysis (absrel)**: Adaptive branch-site random effects likelihood tests
- **Transcriptome Analysis**: Gene expression data for nectar-feeding birds

## Repository Structure

```
nectargenomics/
├── CNEEanalysis_2021/       # Conserved non-coding element analysis
│   ├── 00_inputs/           # Input data and species lists
│   ├── 01a_phylofit/        # Phylogenetic model fitting
│   ├── 01b_cnee_prep/       # CNEE identification and preparation
│   ├── 02_cnee_liftover/    # Coordinate conversion between assemblies
│   └── 03_phyloacc/         # PhyloAcc analysis for accelerated evolution
├── csubst/                  # Convergent substitution analysis
│   ├── 00_inputs/           # Input sequences and species lists
│   └── 01_run_csubst/       # csubst pipeline execution
├── absrel/                  # Positive selection analysis (submodule)
├── Transcriptome_for_NectarGenomics/  # Transcriptome data (submodule)
└── LICENSE                  # GNU GPL v3 license
```

## Species Data

The analyses include data from 45+ bird species, including:
- Multiple hummingbird species (Trochilidae)
- Sunbirds and relatives (Nectariniidae)
- Honeyeaters (Meliphagidae)
- Outgroup species (e.g., chicken - *Gallus gallus*)

Species codes and details can be found in:
- `CNEEanalysis_2021/00_inputs/sp_list.txt`
- `csubst/00_inputs/new_sp_list.txt`

## Main Analyses

### 1. CNEE Analysis

Identifies and analyzes conserved non-coding elements that may be under accelerated evolution in nectar-feeding lineages.

**Pipeline steps:**
1. **PhyloFit**: Fit phylogenetic models to conserved regions
2. **CNEE Preparation**: Identify and annotate CNEEs
3. **LiftOver**: Map CNEEs across genome assemblies
4. **PhyloAcc**: Test for accelerated evolution in target lineages

**Key scripts:**
- `01a_phylofit/03_run_phylofit/` - Phylogenetic model fitting
- `01b_cnee_prep/` - CNEE identification scripts
- `03_phyloacc/` - PhyloAcc analysis for nectar and control lineages

### 2. Convergent Substitution Analysis

Identifies amino acid substitutions that have occurred independently in multiple nectar-feeding lineages.

**Requirements:**
- Python 3.9
- csubst v1.4.0
- IQ-TREE, MAFFT, seqkit
- Nextflow for pipeline execution

**Setup:**
See `csubst/01_run_csubst/README.md` for detailed installation and execution instructions.

### 3. Positive Selection Analysis (absrel)

Tests for positive selection on specific branches using adaptive branch-site random effects likelihood (aBSREL) methods.

This analysis is contained in a separate submodule. See the [absrel repository](https://github.com/osipovarev/absrel) for details.

### 4. Transcriptome Analysis

Gene expression data and analysis for nectar-feeding bird species.

This analysis is contained in a separate submodule. See the [Transcriptome_for_NectarGenomics repository](https://github.com/osipovarev/Transcriptome_for_NectarGenomics) for details.

## Requirements

### Software Dependencies

- **Python**: 3.9+
- **R**: 3.6+ (for statistical analyses and visualization)
- **Conda/Anaconda**: For environment management
- **PhyloAcc**: 2.2.0+ (CNEE analysis)
- **csubst**: 1.4.0 (convergent substitution analysis)
- **IQ-TREE**: 2.2.6+ (phylogenetic inference)
- **MAFFT**: 7.520+ (multiple sequence alignment)
- **seqkit**: (sequence manipulation)
- **Nextflow**: 23.10.0+ (workflow management)
- **SLURM**: For HPC job scheduling (optional)

### Python Packages

```bash
# Core dependencies
numpy
cython
ete3/ete4
pyvolve
matplotlib
biopython
pypdb
pymol-open-source
```

## Installation

### 1. Clone the Repository

```bash
git clone --recurse-submodules https://github.com/maggieMCKO/nectargenomics.git
cd nectargenomics
```

### 2. Set Up Conda Environment

For CNEE analysis (PhyloAcc):
```bash
conda create -n phyloacc220
conda activate phyloacc220
conda install phyloacc
phyloacc.py --version  # Verify PhyloAcc version 2.2.0
```

For convergent substitution analysis:
```bash
conda create --name env_csubst39_14 python=3.9
conda activate env_csubst39_14
conda install iqtree
conda install -c conda-forge pymol-open-source
conda install bioconda::seqkit
conda install bioconda::mafft
pip install git+https://github.com/kfuku52/csubst.git@v1.4.0
pip install https://github.com/etetoolkit/ete/archive/ete4.zip
```

## Usage

### Running CNEE Analysis

1. Prepare input files (MAF alignments, GFF annotations, species tree)
2. Run phylofit to estimate neutral evolution models
3. Identify CNEEs using PhyloP
4. Run PhyloAcc to test for accelerated evolution

See `CNEEanalysis_2021/03_phyloacc/*/README.md` for specific commands.

### Running Convergent Substitution Analysis

```bash
cd csubst/01_run_csubst/nextflow
module load nextflow/23.10.0
conda activate env_csubst39_14
nextflow csubst_mk.nf -profile gwdg -c nextflow.config
```

See `csubst/01_run_csubst/README.md` for detailed instructions.

## Output Files

- **CNEE Results**: Lists of accelerated CNEEs, associated genes, GO enrichment results
- **csubst Results**: Convergent substitution sites, protein structure mappings
- **Phylogenetic Models**: Neutral evolution parameters for each chromosome
- **Statistical Summaries**: P-values, effect sizes, enrichment statistics

## Contributing

This is a research project repository. For questions or collaboration inquiries, please open an issue or contact the repository maintainer.

## Citation

If you use this code or data in your research, please cite:

*(Citation information to be added upon publication)*

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- PhyloAcc developers for the accelerated evolution detection tool
- csubst developers for the convergent substitution analysis framework
- All contributors to the bird genome sequencing projects

## Contact

For questions or issues, please open an issue on GitHub or contact the repository maintainer.

---

**Note**: Some analyses require access to high-performance computing (HPC) resources and large genomic datasets. Execution times and memory requirements vary depending on the number of species and genomic regions analyzed.