#!/usr/bin/bash

#SBATCH -t 0-00:30:00
#SBATCH --mem 40g
#SBATCH -p medium
#SBATCH -n 5
#SBATCH -N 1
#SBATCH --array=1-100

# source ~/default_env.sh
module load conda
source activate R4_gwdg
DATAPATH=/home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/2_GO/1_enrichment
Rscript run_enrichment_perms_gwdg.R galgal6 ${SLURM_ARRAY_TASK_ID} ${SLURM_NTASKS} 50 $DATAPATH
# arg 1 annotation file name,
# arg 2 permutation index ID for slurm batch processing
# arg 3 number of cores
# arg 4 number of permutations
# arg 5 is data path
