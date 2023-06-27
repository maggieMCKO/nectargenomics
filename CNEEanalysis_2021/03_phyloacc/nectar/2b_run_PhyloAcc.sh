#!/bin/bash
#SBATCH -p medium
#SBATCH -n 12             # each task
#SBATCH -N 1              # each task
#SBATCH -t 0-9:00:00      # each task
#SBATCH --mem-per-cpu 6g  # each task
#SBATCH --array 0-363 # ls /scratch2/mko/phyloacc/Target_NeFr40_batches/ | grep -c "batch" then -1
#SBATCH -C "scratch2"
#SBATCH --output=log2_Target_NeFr40/log2_NeFr40_%A_%a.out
#SBATCH --error=log2_Target_NeFr40/log2_NeFr40_%A_%a.err

source "0_setupPhyloAcc.sh"

RUN_num="Target_NeFr40"
PhyloAcc "${scratch_DIR}${RUN_num}_param/run${SLURM_ARRAY_TASK_ID}"
