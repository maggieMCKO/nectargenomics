#!/bin/bash

#### Building a set of consensus CNEEs from literature ####
# based on: https://github.com/sjswuitchik/duck_comp_gen/blob/master/03a_cnee_analysis/01_assemble_cnees.sh

### SETUP
export cneeAna_DIR="../"
export chicken_anno="${cneeAna_DIR}00_inputs/GCF_000002315.6_GRCg6a_genomic.gff" # ori

export replace_chrs_perl="${cneeAna_DIR}01a_phylofit/01_gff_prep/replace_chrs.pl"
export trans_matrix_for_maf="${cneeAna_DIR}01a_phylofit/01_gff_prep/galgal6a_2315.6_acc_chr_matrix.tsv" # (chr)

export merged_bed="bed_outputs/galGal6_all_filtered_merged_5bp.bed" # input
export merged_bed_filtered="bed_outputs/galGal6_all_filtered_merged_5bp_filtered.bed" # output
export merged_bed_filtered_len="cnee_len_ana/galGal6_all_filtered_merged_5bp_filtered.lengths" # output

module purge
module load bedtools/2.29.1

### RUN
# get chicken exons
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.gff.gz
# gunzip GCF_000002315.6_GRCg6a_genomic.gff.gz
awk 'BEGIN{OFS="\t";} {if ($3 ~ /exon/) print $1, $4-1, $5}' ${chicken_anno} | bedtools sort -i - | bedtools merge -i - > bed_outputs/galGal6.exon.bed

# remove small (<50bp) CNEEs
awk '{if ($3-$2 >= 50) print $0}' ${merged_bed} > ${merged_bed_filtered}
awk '{print $3-$2}' ${merged_bed_filtered} > ${merged_bed_filtered_len}

# get CNEEs
bedtools intersect -v -a ${merged_bed_filtered} -b bed_outputs/galGal6.exon.bed > bed_outputs/galGal6_final_merged_CNEEs.bed

# add names
awk 'BEGIN {FS="\t"; OFS="\t"} {$4="CNEE"NR; print}' bed_outputs/galGal6_final_merged_CNEEs.bed > bed_outputs/galGal6_final_merged_CNEEs_named.bed

# replace chr
perl ${replace_chrs_perl} ${trans_matrix_for_maf} bed_outputs/galGal6_final_merged_CNEEs_named.bed > bed_outputs/galGal6_final_merged_CNEEs_named_fixchr_justchr.bed
