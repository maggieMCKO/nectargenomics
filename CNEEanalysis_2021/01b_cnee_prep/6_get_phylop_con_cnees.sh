#!/bin/bash

### SETUP
export cneeAna_DIR="../"
export phyloP_con=phyloP_conserved_cnees_to_keep.tsv

export features_bed=bed_outputs/galGal6_final_merged_CNEEs_named_fixchr_justchr.bed

export filtered_final=bed_outputs/galGal6_final_conserved_CNEEs.bed # final output

### RUN
cut -f4 ${features_bed}  | grep -Fx -f ${phyloP_con}> ${filtered_final}
wc -l ${filtered_final}
