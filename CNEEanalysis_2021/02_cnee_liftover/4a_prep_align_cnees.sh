#!/bin/bash
#SBATCH --time=00:30:00 
#SBATCH -c 1
#SBATCH --mem=5g
#SBATCH --partition=fat
#SBATCH -o log4a_prepalign_%A.out
#SBATCH -e log4a_prepalign_%A.err
#SBATCH -J prepalign

# WD: /home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/02_cnee_liftover
source "0_setup.sh"
mkdir -p log5d

export dir_cneeconcatfa=liftover_proc/4_concatfa/

all_cnees_tab_from="${dir_cneeconcatfa}all_cnees.tab"
all_cnees_tab_to="${scratch_DIR}all_cnees.tab"

# use scratch system
cp ${all_cnees_tab_from} ${all_cnees_tab_to}
echo 'copy tab done'

## 1. aligned (galGal6)
cut -f4,4 ${filtered_final} | split -a 3 -d  -l 1000 - ${aligned_DIR}batch
echo 'split align done'
