#!/bin/bash
#SBATCH --time=04:10:00 # 211mins
#SBATCH -c 1
#SBATCH --mem=4g
#SBATCH --partition=medium
#SBATCH -o log3_generateCNEEfa_%A.out
#SBATCH -e log3_generateCNEEfa_%A.err
#SBATCH -J generateCNEEfa
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mcko@orn.mpg.de
#SBATCH -C "scratch|scratch2"

# adapted from https://github.com/sjswuitchik/duck_comp_gen/blob/2bf589d198122fd9ae69337a30e62ae303c7c9cb/03a_cnee_analysis/02_align_cnees.sh

####### SETUP
source "0_setup.sh"


# use scratch: setup in 0_setup.sh


all_cnees_tab="all_cnees.tab"
all_cnees_tab_from="${tmp_pslMappedTogalGal_dir}${all_cnees_tab}"
all_cnees_tab_to="${MYSCRATCH}${all_cnees_tab}"

####### RUN
## 0. use scratch system
cd ${MYSCRATCH}
cp ${all_cnees_tab_from} ${all_cnees_tab_to}

## 1. aligned
cut -f4,4 ${tmp_pslMappedTogalGal_dir}galGal6_cnees_parsed_liftover.bed | split -a 3 -d  -l 1000 - batch
mv batch* ${aligned_DIR}

## 2. unaligned
awk -v targetDir="$unaligned_DIR" '{print ">"$1 >> targetDir$2".fa"; print $3 >> targetDir$2".fa"; close(targetDir$2".fa")}' ${all_cnees_tab_to}
