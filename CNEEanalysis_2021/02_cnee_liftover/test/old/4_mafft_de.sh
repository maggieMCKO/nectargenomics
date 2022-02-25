#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 0-02:30 # for each tak; 30 mins per task
#SBATCH --output=log4_%A_%a.out
#SBATCH --error=log4_%A_%a.err
#SBATCH -p medium
#SBATCH --mem=2000
#SBATCH --array=0-316 ###!!!!!!!! output of ls $aligned_DIR | wc -l; then minus 1
#SBATCH -C "scratch|scratch2"

# adapted from https://github.com/sjswuitchik/duck_comp_gen/blob/2bf589d198122fd9ae69337a30e62ae303c7c9cb/03a_cnee_analysis/03_mafft.sh

####### SETUP
source "0_setup.sh"

# use scratch: setup in 0_setup.sh

export trimalAL_PATH="$HOME/Tools/trimal/source/trimal"
module load mafft/7.304
module load anaconda2/2019.10
source activate py27 # for using perl v5

####### RUN
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
