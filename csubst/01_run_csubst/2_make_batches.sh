#!/bin/bash


# wd: /home/mpg08/mko/Nectar/analysis/csubst/01_run_csubst

export wd="$HOME/Nectar/analysis/csubst/01_run_csubst/"
export input_dir="$HOME/Nectar/analysis/csubst/00_inputs/"
export gene_list="${input_dir}isoforms.one2ones.checked.tsv"
export trans_list="${input_dir}checked.transcripts.lst"

export scratch_DIR="/scratch/users/$USER/csubst/"
export batch_DIR="${wd}batches/"

mkdir -p ${scratch_DIR} ${batch_DIR}


num_per_batch=20
split -d -a 4 -l ${num_per_batch} ${trans_list} ${batch_DIR}

