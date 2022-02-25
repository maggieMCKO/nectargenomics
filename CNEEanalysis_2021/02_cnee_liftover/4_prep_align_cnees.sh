#!/bin/bash
#SBATCH --time=04:00:00 # 186-211mins
#SBATCH -c 1
#SBATCH --mem=4g
#SBATCH --partition=medium
#SBATCH -o log4_prepalign_%A.out
#SBATCH -e log4_prepalign_%A.err
#SBATCH -J prepalign
##SBATCH --mail-type=ALL
##SBATCH --mail-user=mcko@orn.mpg.de
#SBATCH -C "scratch2"

# WD: /home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/02_cnee_liftover
source "0_setup.sh"

export dir_cneeconcatfa=liftover_proc/4_concatfa/

all_cnees_tab_from="${dir_cneeconcatfa}all_cnees.tab"
all_cnees_tab_to="${scratch_DIR}all_cnees.tab"

# use scratch system
cp ${all_cnees_tab_from} ${all_cnees_tab_to}

## 1. aligned (galGal6)
cut -f4,4 ${filtered_final} | split -a 3 -d  -l 1000 - ${aligned_DIR}batch

## 2. unaligned
awk -v targetDir="$unaligned_DIR" '{print ">"$1 >> targetDir$2".fa"; print $3 >> targetDir$2".fa"; close(targetDir$2".fa")}' ${all_cnees_tab_to}

mkdir -p log5d
