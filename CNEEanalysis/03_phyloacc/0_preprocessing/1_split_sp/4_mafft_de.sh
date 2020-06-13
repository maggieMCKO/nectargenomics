#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 0-01:30
#SBATCH -p medium
#SBATCH --mem=2000
#SBATCH --array=0-323 ###!!!!!!!! length of wc -l aligned
#SBATCH -C "scratch2"

# adapted from https://github.com/sjswuitchik/duck_comp_gen/blob/2bf589d198122fd9ae69337a30e62ae303c7c9cb/03a_cnee_analysis/03_mafft.sh

####### SETUP
export cneeAna_DIR="/home/mpg08/mko/Nectar/analysis/CNEEanalysis/"
export WD="${cneeAna_DIR}03_phyloacc/preprocessing/1_way2/"
export specieslist="${WD}Sp.list"
export trans_matrix2="${cneeAna_DIR}00_maf_gff_preprocessing/galgal6a_2315.6_acc5.tsv"

# use scratch
MYSCRATCH="/scratch/${USER}/CNEE/" 
aligned_DIR="${MYSCRATCH}aligned/"
unaligned_DIR="${MYSCRATCH}unaligned/"

export trimalAL_PATH="/usr/users/mko/Tools/trimal/source/trimal"
module load MAFFT
module load conda
source activate py27 # for using perl v5

####### RUN
aligned_DIR="${MYSCRATCH}aligned/"
cd ${aligned_DIR}

printf -v BATCH "%03d" $SLURM_ARRAY_TASK_ID
mkdir -p batch${BATCH}_output

for CNEE in $(cat batch$BATCH);
do
  if [ ! -s batch${BATCH}_output/$CNEE.aligned.fa ]
  then
    ginsi ${unaligned_DIR}${CNEE}.fa | perl -p -e 'if (!/^>/) {s/[Nn]/-/g}' > temp_${BATCH}.fa
    ${trimalAL_PATH} -in temp_${BATCH}.fa -noallgaps -out batch${BATCH}_output/$CNEE.aligned.fa -fasta
    rm temp_${BATCH}.fa
  fi
done
