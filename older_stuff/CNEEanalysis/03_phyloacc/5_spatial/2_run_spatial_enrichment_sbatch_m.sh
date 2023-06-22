#!/usr/bin/bash

#SBATCH -t 0-01:00:00
#SBATCH --mem 40g
#SBATCH -p medium
#SBATCH -n 1
#SBATCH -N 1
##SBATCH --array=1-100


# source ~/default_env.sh
module load conda
source activate R402_gwdg
Rscript run_spatial_enrichment_gwdg_m.R 
