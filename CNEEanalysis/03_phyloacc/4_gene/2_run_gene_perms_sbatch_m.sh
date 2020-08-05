#!/usr/bin/bash

#SBATCH -t 0-00:10:00
#SBATCH --mem 40g
#SBATCH -p medium
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --array=1-100


# source ~/default_env.sh
module load conda
source activate R402_gwdg
DATAPATH=/home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/4_gene
Rscript run_gene_perms_gwdg_m.R galgal6 ${SLURM_ARRAY_TASK_ID} ${SLURM_NTASKS} 100 $DATAPATH
# arg 1 annotation file name,
# arg 2 permutation index ID for slurm batch processing
# arg 3 number of cores
# arg 4 number of permutations
# arg 5 is data path
