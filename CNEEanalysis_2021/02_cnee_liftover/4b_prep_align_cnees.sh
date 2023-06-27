#!/bin/bash
#SBATCH --time=12:00:00 # 186-211mins, once 569mins!
#SBATCH -c 2
#SBATCH --mem=8g
#SBATCH --partition=fat
#SBATCH -o log4b_prepalign_%A.out
#SBATCH -e log4b_prepalign_%A.err
#SBATCH -J prepalign

# WD: /home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/02_cnee_liftover
source "0_setup.sh"

export all_cnees_tab_to="${scratch_DIR}all_cnees.tab"

## 2. unaligned
awk -v targetDir="$unaligned_DIR" '{print ">"$1 >> targetDir$2".fa"; print $3 >> targetDir$2".fa"; close(targetDir$2".fa")}' ${all_cnees_tab_to}
