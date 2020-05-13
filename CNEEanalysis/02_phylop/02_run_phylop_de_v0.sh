#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4g
#SBATCH --partition=medium
#SBATCH -o phylop_%A.out
#SBATCH -e phylop_%A.err
#SBATCH -J phylop
#SBATCH --mail-type=ALL        # BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mcko@orn.mpg.de

##### SETUP INPUT and OUTPUT
cneeAna_DIR="/home/mpg08/mko/Nectar/analysis/CNEEanalysis/"

maf_DIR="${cneeAna_DIR}MAF/"
input_maf="${maf_DIR}multi.anno.maf"

phylofit_DIR="${cneeAna_DIR}01_phylofit/01_msa3/"
neutral_model="${phylofit_DIR}nonconserved_4d.mod"

WD="${cneeAna_DIR}02_phylop/"
features_bed="${WD}galGal6_final_merged_CNEEs_named_fixchr2.bed"
features_out="${WD}features.out"

##### SETUP programs
export phastPATH="/home/mpg08/mko/Tools/phast/bin/"
export phyloP_PATH="${phastPATH}phyloP"

##### RUN
### Run phyloP with the features option.
# USAGE: phyloP [OPTIONS] tree.mod [alignment] > out
${phyloP_PATH} --method LRT --mode CONACC --msa-format MAF --features ${features_bed} ${neutral_model} ${input_maf} > ${features_out}
