#!/bin/bash

#### Building a set of consensus CNEEs from literature ####
# based on: https://github.com/sjswuitchik/duck_comp_gen/blob/master/03a_cnee_analysis/01_assemble_cnees.sh

### SETUP
export cneeAna_DIR="$HOME/Nectar/analysis/CNEEanalysis_2021/"
export chicken_anno="${cneeAna_DIR}00_inputs/GCF_000002315.6_GRCg6a_genomic.gff" # ori

export replace_chrs_perl="${cneeAna_DIR}01a_phylofit/01_gff_prep/replace_chrs.pl"
export trans_matrix_for_maf="${cneeAna_DIR}01a_phylofit/01_gff_prep/galgal6a_2315.6_acc_chr_matrix.tsv" # (chr)

### OUTPUT
export features_bed="${cneeAna_DIR}01b_cnee_prep/galGal6_final_merged_CNEEs_named_fixchr_justchr.bed"

module purge
module load bedtools/2.29.1

### RUN
# get ce lengths
awk '{print $3-$2, "\tLowe"}' "${cneeAna_DIR}00_inputs/cnee_bedfiles/galGal6_Lowe_merged.bed" >> ce.lengths
awk '{print $3-$2, "\tSackton"}' "${cneeAna_DIR}00_inputs/cnee_bedfiles/galGal6_Sackton_merged.bed"  >> ce.lengths
awk '{print $3-$2, "\tUCSC"}' "${cneeAna_DIR}00_inputs/cnee_bedfiles/galGal6_UCSC_merged.bed" >> ce.lengths
awk '{print $3-$2, "\tCraig"}' "${cneeAna_DIR}00_inputs/cnee_bedfiles/galGal6_Craig_merged.bed" >> ce.lengths

# merge cnee sets
for d in 0 1 2 3 4 5 6 7 8 9 10
do
cat ${cneeAna_DIR}00_inputs/cnee_bedfiles/galGal6_*_merged.bed | bedtools sort -i - | bedtools merge -i - -d $d | awk '{if ($3-$2 >= 50) print $0}' > galGal6_all_merged_${d}bp.bed
awk '{print $3-$2}' galGal6_all_merged_${d}bp.bed > galGal6_all_merged_${d}bp.lengths
done
