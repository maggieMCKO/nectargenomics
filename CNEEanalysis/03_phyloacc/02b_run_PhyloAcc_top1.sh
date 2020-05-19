#!/bin/bash
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 2-00:00
#SBATCH --mem 80G # 32000 not enough
#SBATCH --array 0-1877

source /n/holylfs/LABS/informatics/nectar_genomes/analyses/CNEEanalysis/03_phyloacc/setupPhyloAcc.sh
PhyloAcc top1_param/run${SLURM_ARRAY_TASK_ID}
