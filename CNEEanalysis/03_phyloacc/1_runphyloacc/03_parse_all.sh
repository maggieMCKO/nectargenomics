#!/bin/bash
# based on https://github.com/sjswuitchik/duck_comp_gen/blob/master/03a_cnee_analysis/06_PhyloAcc/parse_all.sh
## Concat output from PhyloAcc ##

source /home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/1_run_phyloacc/setupPhyloAcc.sh
cd ${MYSCRATCH}
cd top1_outs

VERSION=1
concat_output_PATH="${WD}concat_output.sh"


${concat_output_PATH} elem_lik
mv elem_lik_combined.txt ../elem_lik_combined_top1_${VERSION}.txt
cp ../elem_lik_combined_top1_${VERSION}.txt ${WD}

${concat_output_PATH} rate_postZ_M2
mv rate_postZ_M2_combined.txt ../rate_postZ_M2_combined_top1_${VERSION}.txt
cp ../rate_postZ_M2_combined_top1_${VERSION}.txt ${WD}
