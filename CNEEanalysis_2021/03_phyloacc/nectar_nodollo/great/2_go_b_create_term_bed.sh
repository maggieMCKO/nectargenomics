#!/bin/bash
#SBATCH --time=0:30:00
#SBATCH --ntasks=24
#SBATCH --mem-per-cpu=5g
#SBATCH --partition=medium
#SBATCH -o log2_create_go_term_bed_%A.o
#SBATCH -e log2_create_go_term_bed_%A.e
#SBATCH -J create_go_term


# wd: /home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/05_perm_phyloacc/1_rmPar/1_great


export great_dir="/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/03_phyloacc/2024ver/nectar_nodollo/great/"
export regDoms1="${great_dir}galGal6_gene_tss_great_basalPlusExtension_proteincoding.bed"
# export regDoms2="${great_dir}galGal6_gene_tss_great_oneClosest_proteincoding.bed"
# export regDoms3="${great_dir}galGal6_gene_tss_great_twoClosest_proteincoding.bed"


# I. create bed for each GO term
module load rev/23.12 anaconda3/2023.09-0
source activate R432 

Rscript 2_go_a_subset.R 24 ${regDoms1}
# 1: n core; 2: full path of great output

# Rscript 2_go_a_subset.R 24 ${regDoms2}
# Rscript 2_go_a_subset.R 24 ${regDoms3}
