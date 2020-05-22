#!/bin/bash
#SBATCH --time=0-01:30:00
#SBATCH -n 1 # or --ntasks=1
#SBATCH --mem=20g
#SBATCH --partition=medium
#SBATCH -o phylopBC_%A.out
#SBATCH -e phylopBC_%A.err
#SBATCH -J phylopBC
#SBATCH --mail-type=ALL        # BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mcko@orn.mpg.de

source setup_paths.sh
# WD="${cneeAna_DIR}02_phylop/subtree/"

tmp_branch=$1
echo -e "tmp_branch: ${tmp_branch}"

tmp_bed=$2
echo -e "tmp_bed: ${tmp_bed}"

chr=$3
echo -e "tmp_chr: ${chr}"

tmp_maf=$4
echo -e "tmp_maf: ${tmp_maf}"

method=LRT
echo -e "method: ${method}"

mode=CON
echo -e "mode: ${mode}"

mkdir -p out_branch
features_out="out_branch/${chr}_${method}_${mode}_${tmp_branch}_cnee.wig"
echo -e "features_out: ${features_out}" # in WD

##### RUN
# USAGE: phyloP [OPTIONS] tree.mod [alignment] > out
${phyloP_PATH} --method ${method} --branch ${tmp_branch} --mode ${mode} --msa-format MAF --features ${tmp_bed} ${neutral_model_named} ${tmp_maf} > ${features_out}
