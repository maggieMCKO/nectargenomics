#!/usr/bin/bash

#SBATCH -t 0-00:30:00
#SBATCH --mem 40g
#SBATCH -p medium
#SBATCH -n 1
#SBATCH -N 1

# source ~/default_env.sh
module load conda
source activate R402_gwdg

Rscript run_go_ecdf_calculation_gwdg.R
