#!/bin/bash

export Ref_species="galGal6"
export specieslist="../00_inputs/sp_list.txt"
export filtered_final="../01b_cnee_prep/bed_outputs/galGal6_final_conserved_CNEEs.bed"

export scratch_DIR="/scratch/users/$USER/phyloacc/"
export aligned_DIR="${scratch_DIR}aligned/"
export unaligned_DIR="${scratch_DIR}unaligned/"
mkdir -p ${scratch_DIR} ${aligned_DIR} ${unaligned_DIR}
