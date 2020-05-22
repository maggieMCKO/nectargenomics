#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --mem=20g
#SBATCH --partition=medium
#SBATCH -o phylop_%A.out
#SBATCH -e phylop_%A.err
#SBATCH -J phylop
#SBATCH --mail-type=ALL        # BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mcko@orn.mpg.de

source setup_paths.sh


tmp_bed=$1
echo -e "tmp_bed: ${tmp_bed}" # full path

chr=$2
echo -e "tmp_chr: ${chr}"

tmp_maf=$3
echo -e "tmp_maf: ${tmp_maf}"

method=LRT
echo -e "method: ${method}"

mode=CON
echo -e "mode: ${mode}"

mkdir -p out

features_out="out/${chr}_${method}_${mode}_cnee.wig"
echo -e "features_out: ${features_out}" # current WD

##### RUN
# USAGE: phyloP [OPTIONS] tree.mod [alignment] > out
${phyloP_PATH} --method ${method} --mode ${mode} --msa-format MAF --features ${tmp_bed} ${neutral_model} ${tmp_maf} > ${features_out}
