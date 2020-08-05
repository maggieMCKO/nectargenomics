#!/usr/bin/bash

#SBATCH -t 0-00:30:00
#SBATCH --mem 160g
#SBATCH -p fat
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -o ecdf_cal_%A.o
#SBATCH -e ecdf_cal_%A.e

# source ~/default_env.sh
module load conda
source activate R402_gwdg

Rscript run_gene_ecdf_calculation_gwdg_m.R
