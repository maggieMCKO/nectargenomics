#!/bin/bash
#SBATCH --time=01:30:00
#SBATCH --ntasks=1
#SBATCH --mem=20g
#SBATCH --partition=medium
#SBATCH -o phylop_logs/phylop_%A.out
#SBATCH -e phylop_logs/phylop_%A.err
#SBATCH -J phylop

export neutral_model="../01a_phylofit/03_run_phylofit/nonconserved_4d.mod"

export phastPATH="/home/mpg08/mko/Tools/phast/bin/"
export phyloP_PATH="${phastPATH}phyloP"


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

tmp_out_dir="phyloP_outputs/"
mkdir -p ${tmp_out_dir}
features_out="${tmp_out_dir}${chr}_${method}_${mode}_cnee.tsv"
echo -e "features_out: ${features_out}" # current WD

##### RUN
### Run phyloP with the features option.
# USAGE: phyloP [OPTIONS] tree.mod [alignment] > out
${phyloP_PATH} --method ${method} --mode ${mode} --msa-format MAF --features ${tmp_bed} ${neutral_model} ${tmp_maf} > ${features_out}
