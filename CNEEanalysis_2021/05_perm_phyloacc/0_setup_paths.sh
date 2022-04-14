#!/bin/bash
##### SETUP INPUT and OUTPUT
export cneeAna_DIR="$HOME/Nectar/analysis/CNEEanalysis_2021/"
export input_DIR="${cneeAna_DIR}00_inputs/"
export input_tree="${input_DIR}NectarTree02.nw"

# export chicken_anno="${input_DIR}GCF_000002315.6_GRCg6a_genomic_chrtabFixed2.gff" # (galgal.chrx) # seems dont work
export chicken_anno="${input_DIR}GCF_000002315.6_GRCg6a_genomic_chrtabFixed_justchr.gff" # just chr

export features_bed="${cneeAna_DIR}01b_cnee_prep/bed_outputs/galGal6_final_conserved_CNEEs.bed"
