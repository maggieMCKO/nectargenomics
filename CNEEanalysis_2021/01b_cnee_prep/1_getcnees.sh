#!/bin/bash

#### Building a set of consensus CNEEs from literature ####
# based on: https://github.com/sjswuitchik/duck_comp_gen/blob/master/03a_cnee_analysis/01_assemble_cnees.sh

### SETUP
# export cneeAna_DIR="$HOME/Nectar/analysis/CNEEanalysis_2021/"
export cneeAna_DIR="../"

module purge
module load bedtools/2.29.1

### RUN
mkdir -p cnee_len_ana
mkdir -p bed_outputs

awk '{print $3-$2, "\tLowe"}' "${cneeAna_DIR}00_inputs/cnee_bedfiles/galGal6_Lowe_merged.bed" >> cnee_len_ana/ce.lengths
awk '{print $3-$2, "\tSackton"}' "${cneeAna_DIR}00_inputs/cnee_bedfiles/galGal6_Sackton_merged.bed"  >> cnee_len_ana/ce.lengths
awk '{print $3-$2, "\tUCSC"}' "${cneeAna_DIR}00_inputs/cnee_bedfiles/galGal6_UCSC_merged.bed" >> cnee_len_ana/ce.lengths
awk '{print $3-$2, "\tCraig"}' "${cneeAna_DIR}00_inputs/cnee_bedfiles/galGal6_Craig_merged.bed" >> cnee_len_ana/ce.lengths

# remove cnees larger than 5000 bp, then merge cnee sets
awk '{if ($3-$2 <= 5000) print $0}' "${cneeAna_DIR}00_inputs/cnee_bedfiles/galGal6_Craig_merged.bed" > bed_outputs/galGal6_Craig_filtered.bed
awk '{if ($3-$2 <= 5000) print $0}' "${cneeAna_DIR}00_inputs/cnee_bedfiles/galGal6_Lowe_merged.bed" > bed_outputs/galGal6_Lowe_filtered.bed
awk '{if ($3-$2 <= 5000) print $0}' "${cneeAna_DIR}00_inputs/cnee_bedfiles/galGal6_Sackton_merged.bed"  > bed_outputs/galGal6_Sackton_filtered.bed
awk '{if ($3-$2 <= 5000) print $0}' "${cneeAna_DIR}00_inputs/cnee_bedfiles/galGal6_UCSC_merged.bed" > bed_outputs/galGal6_UCSC_filtered.bed

# merge cnee sets
for d in 0 1 2 3 4 5 6 7 8 9 10
do
  out_bed="bed_outputs/galGal6_all_filtered_merged_${d}bp.bed"
  out_len="cnee_len_ana/galGal6_all_filtered_merged_${d}bp.lengths"
  cat bed_outputs/galGal6_*_filtered.bed | bedtools sort -i - | bedtools merge -i - -d $d > ${out_bed}
  awk '{print $3-$2}' ${out_bed} > ${out_len}
done
