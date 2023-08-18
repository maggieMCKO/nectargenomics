#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --ntasks=24
#SBATCH --mem=160g
#SBATCH --partition=fat
#SBATCH -o log_concat1_%A.o
#SBATCH -e log_concat1_%A.e
#SBATCH -J concat1

# wd: /home/mpg08/mko/Nectar/analysis/pcoc/00_inputs

module load anaconda2/2019.10
source activate R402_2022_2

Rscript 5a_concatenate_results_v2.R /home/mpg08/mko/Nectar/analysis/csubst/01_run_csubst/out_v3/ run3 24
