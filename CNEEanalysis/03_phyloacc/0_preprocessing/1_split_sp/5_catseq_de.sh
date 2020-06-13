#!/bin/bash
#SBATCH --time=01:10:00
#SBATCH -c 1
#SBATCH --mem=4g
#SBATCH --partition=medium
#SBATCH -o log5_concatseq_%A.out
#SBATCH -e log5_concatseq_%A.err
#SBATCH -J concatseq
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mcko@orn.mpg.de
#SBATCH -C "scratch2"

# adapted from https://github.com/sjswuitchik/duck_comp_gen/blob/2bf589d198122fd9ae69337a30e62ae303c7c9cb/03a_cnee_analysis/04_catseq.sh

## Concatenation for CNEE analysis ##
## https://github.com/ChrisCreevey/catsequences ##
## NB: run from git clone  --branch seqname https://github.com/harvardinformatics/catsequences.git until PR is fulfilled ##

####### SETUP
export cneeAna_DIR="/home/mpg08/mko/Nectar/analysis/CNEEanalysis/"
export WD="${cneeAna_DIR}03_phyloacc/preprocessing/1_way2/"
export specieslist="${WD}Sp.list"
export tmp_pslMappedTogalGal_dir="${WD}CneeMappedToGg6_PSL_DIR/"

# use scratch
MYSCRATCH="/scratch/${USER}/CNEE/" 
aligned_DIR="${MYSCRATCH}aligned/"
unaligned_DIR="${MYSCRATCH}unaligned/"

export catsequences_PATH="/usr/users/mko/Tools/catsequences/catsequences"

####### RUN
aligned_DIR="${MYSCRATCH}aligned/"
cd ${aligned_DIR}

# concat
find . -name '*.fa' > list
${catsequences_PATH} list

mv allseqs.fas galloseq.fa
mv allseqs.partitions.txt galloseq.partitions.txt

# clean up ?s
cat galloseq.fa | perl -p -e 's/[?]/-/g' > galloseq_gapFixed.fa

# need to make part.txt into a bed with CNEE-start-end
awk 'BEGIN{FS="="; OFS="\t"} {split($1,a,"_output/"); print a[2],$2}' galloseq.partitions.txt | sed 's/;$//' | sed 's/\.aligned\.fa//g' | tr -s '\t'| awk 'BEGIN{FS="\t"; OFS="\t"} {split($2,a,"-"); print $1,a[1],a[2]}' | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $3}' > galloseq.part.bed

# copy to WD
cp galloseq.partitions.txt ${tmp_pslMappedTogalGal_dir}
cp galloseq_gapFixed.fa ${tmp_pslMappedTogalGal_dir}
cp galloseq.part.bed ${tmp_pslMappedTogalGal_dir}
