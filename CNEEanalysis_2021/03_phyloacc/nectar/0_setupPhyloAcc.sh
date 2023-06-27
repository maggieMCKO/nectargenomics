#!/bin/bash

module purge
module load gcc/9.3.0 gsl/2.5 phyloacc/1.0
export interface_PATH="$HOME/Tools/PhyloAcc-interface/src/interface/phyloacc_interface.py"
export PATH=/usr/product/applsw/phyloacc-1.0/install/bin:${PATH} # PhyloAcc
export PATH=/usr/product/applsw/phyloacc-1.0/install.v2_gbgc/bin/:${PATH} # PhyloAcc_gBGC
# because gwdg didn't configure well

# python ${interface_PATH} -h

export scratch_DIR="/scratch/users/$USER/phyloacc/"
export aligned_DIR="${scratch_DIR}aligned/"
export unaligned_DIR="${scratch_DIR}unaligned/"

export neutral_model_named="../../01a_phylofit/03_run_phylofit/nonconserved_4d_named.mod"
export input_cnee_fa="../../02_cnee_liftover/liftover_proc/5_aligned_cnees/galloseq_gapFixed.fa"
export input_part_bed="../../02_cnee_liftover/liftover_proc/5_aligned_cnees/galloseq.part.bed"

export MYSCRATCH_input="${scratch_DIR}input/"
export MYSCRATCH_output="${scratch_DIR}output/"
mkdir -p ${MYSCRATCH_input} ${MYSCRATCH_output}
export scratch_neutral_model_named="${MYSCRATCH_input}nonconserved_4d_named.mod"
export scratch_input_cnee_fa="${MYSCRATCH_input}galloseq_gapFixed.fa"
export scratch_input_part_bed="${MYSCRATCH_input}galloseq.part.bed"
