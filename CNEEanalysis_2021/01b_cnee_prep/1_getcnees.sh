#!/bin/bash

#### Building a set of consensus CNEEs from literature ####
# based on: https://github.com/sjswuitchik/duck_comp_gen/blob/master/03a_cnee_analysis/01_assemble_cnees.sh

### SETUP
export cneeAna_DIR="$HOME/Nectar/analysis/CNEEanalysis_2021/"

module purge
module load bedtools/2.29.1

### RUN
mkdir -p cnee_len_ana
mkdir -p bed_outputs

get cnee_len_ana/ce.lengths
awk '{print $3-$2, "\tLowe"}' "${cneeAna_DIR}00_inputs/cnee_bedfiles/galGal6_Lowe_merged.bed" >> cnee_len_ana/ce.lengths
awk '{print $3-$2, "\tSackton"}' "${cneeAna_DIR}00_inputs/cnee_bedfiles/galGal6_Sackton_merged.bed"  >> cnee_len_ana/ce.lengths
awk '{print $3-$2, "\tUCSC"}' "${cneeAna_DIR}00_inputs/cnee_bedfiles/galGal6_UCSC_merged.bed" >> cnee_len_ana/ce.lengths
awk '{print $3-$2, "\tCraig"}' "${cneeAna_DIR}00_inputs/cnee_bedfiles/galGal6_Craig_merged.bed" >> cnee_len_ana/ce.lengths

# merge cnee sets
for d in 0 1 2 3 4 5 6 7 8 9 10
do
cat ${cneeAna_DIR}00_inputs/cnee_bedfiles/galGal6_*_merged.bed | bedtools sort -i - | bedtools merge -i - -d $d > bed_outputs/galGal6_all_merged_${d}bp.bed

awk '{print $3-$2}' bed_outputs/galGal6_all_merged_${d}bp.bed > cnee_len_ana/galGal6_all_merged_${d}bp.lengths
done

# remove cnees larger than 1000 then merge cnee sets
awk '{if ($3-$2 <= 1000) print $0}' "${cneeAna_DIR}00_inputs/cnee_bedfiles/galGal6_Lowe_merged.bed" > bed_outputs/galGal6_Lowe_merged_eqless1000.bed
awk '{if ($3-$2 <= 1000) print $0}' "${cneeAna_DIR}00_inputs/cnee_bedfiles/galGal6_Sackton_merged.bed"  > bed_outputs/galGal6_Sackton_merged_eqless1000.bed
awk '{if ($3-$2 <= 1000) print $0}' "${cneeAna_DIR}00_inputs/cnee_bedfiles/galGal6_UCSC_merged.bed" > bed_outputs/galGal6_UCSC_merged_eqless1000.bed
awk '{if ($3-$2 <= 1000) print $0}' "${cneeAna_DIR}00_inputs/cnee_bedfiles/galGal6_Craig_merged.bed" > bed_outputs/galGal6_Craig_merged_eqless1000.bed

for d in 0 1 2 3 4 5 6 7 8 9 10
do
cat bed_outputs/galGal6_*_merged_eqless1000.bed | bedtools sort -i - | bedtools merge -i - -d $d > bed_outputs/galGal6_all_merged_eqless1000_${d}bp.bed

awk '{print $3-$2}' bed_outputs/galGal6_all_merged_eqless1000_${d}bp.bed > cnee_len_ana/galGal6_all_merged_eqless1000_${d}bp.lengths
done

# | awk '{if ($3-$2 >= 50) print $0}'
