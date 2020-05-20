#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH --partition=medium
#SBATCH -o MAFtoFA_%A.out
#SBATCH -e MAFtoFA_%A.err
#SBATCH -J MAFtoFA
#SBATCH --mail-type=ALL        # BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mcko@orn.mpg.de

##### SETUP INPUT and OUTPUT
cneeAna_DIR="/home/mpg08/mko/Nectar/analysis/CNEEanalysis/"

maf_DIR="${cneeAna_DIR}MAF/"
input_maf="${maf_DIR}multi.anno.maf"

FA_out_DIR="${maf_DIR}fasta/"


##### SETUP programs
export phastPATH="/home/mpg08/mko/Tools/phast/bin/"
export msa_view_PATH="${phastPATH}msa_view"

##### RUN
Stem=$(basename $input_maf)
Stem_FA=${Stem/.maf/.fasta}
FA_out="${FA_out_DIR}${Stem_FA}"

# USAGE: maf_parse [OPTIONS] <infile>
${msa_view_PATH} --in-format MAF --out-format FASTA ${input_maf} > ${FA_out}
