#!/bin/bash
#SBATCH -p medium
#SBATCH -n 12
#SBATCH -N 1
#SBATCH -t 0-10:00 # time for a single job
#SBATCH --mem 8G # used 4.72
#SBATCH --array 0-161
#SBATCH -C "scratch2"

source /home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/1_run_phyloacc/setupPhyloAcc.sh

PhyloAcc "${MYSCRATCH}top1_param/run${SLURM_ARRAY_TASK_ID}"
