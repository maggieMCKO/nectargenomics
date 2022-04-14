#!/bin/bash


module purge
module load bedtools/2.29.1

old=/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021_old/03_phyloacc/0_preprocessing/1_prepare_cnees/galGal6_final_merged_CNEEs_named_fixchr_justchr.bed
cleaned=/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/01b_cnee_prep/bed_outputs/galGal6_final_conserved_CNEEs.bed

bedtools intersect -loj \
    -a ${old} \
    -b ${cleaned} \
    -names cleaned \
    -sorted > 'comparing_cleaned_cnees.bed'