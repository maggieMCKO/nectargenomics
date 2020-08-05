#!/usr/bin/bash

#SBATCH -t 0-00:15:00
#SBATCH --mem 20g
#SBATCH -p medium
#SBATCH -n 5
#SBATCH -N 1
#SBATCH --array=1-100

# source ~/default_env.sh
module load conda
source activate R402_gwdg
DATAPATH=/home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/3_convergencePerm
Rscript run_convergence_perms_gwdg_m.R ${SLURM_ARRAY_TASK_ID} ${SLURM_NTASKS} 50 $DATAPATH

# arg 1 permutation index ID for slurm batch processing
# arg 3 number of cores
# arg 4 number of permutations
# arg 5 is data path
