#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --ntasks=4
#SBATCH --mem=20g
#SBATCH --partition=fat
#SBATCH -o log_concat_rnaseq_%A.o
#SBATCH -e log_concat_rnaseq_%A.e
#SBATCH -J concat_rnaseq

# wd: /home/mpg08/mko/Nectar/analysis/rnaseq/00_inputs

module load anaconda2/2019.10
source activate R402_2022

Rscript 1a_concatenate_rnaseq_results.R /home/mpg08/mko/Nectar/analysis/rnaseq/00_inputs/Kallisto_quntification_all_lib
