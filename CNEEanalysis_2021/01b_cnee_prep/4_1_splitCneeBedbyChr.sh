#!/bin/bash

export cneeAna_DIR="../"
export trans_matrix_for_maf="${cneeAna_DIR}01a_phylofit/01_gff_prep/galgal6a_2315.6_acc_chr_matrix.tsv" # (chr)

export features_bed=bed_outputs/galGal6_final_merged_CNEEs_named_fixchr_justchr.bed

# 1. split cnee bed by chromosome
mkdir -p cnee_byChr
mkdir -p phylop_logs
splitbedbyChr (){
    Chr=$1
    echo $Chr

    # input bed
    tmp_bed="cnee_byChr/galGal6_${Chr}_justchr.bed"

    awk '/'$Chr'\t/' ${features_bed} > ${tmp_bed}
    # line_bed=$(wc -l $tmp_bed)
    # echo -e "bed: ${line_bed}"
    # line_gff_cut=$(echo ${line_bed} | cut -d ' ' -f1 )
}

export -f splitbedbyChr
cut -f2 ${trans_matrix_for_maf} | parallel splitbedbyChr {}

# bash 4_1_splitCneeBedbyChr.sh &> splitcnee.log
# to check: ls cnee_byChr | xargs wc -l
