# https://github.com/sjswuitchik/duck_comp_gen/blob/master/03a_cnee_analysis/06_PhyloAcc/parse_all.sh
## Concat output from PhyloAcc ##

### setup
source "0_setupPhyloAcc.sh"
export WD=$(pwd)
export concat_output_PATH="${WD}/0_concat_output.sh" # need absolute path

# WD: /home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/03_phyloacc/oldway

# Runs=("Target_Ne40" "Target_NeFr50" "Target_NeFr50cli" "Control40" "Control50")
# Runs=("Target_NeFr40_aca" "Target_NeFr40_cliaca" "Target_NeFr40_songbirdFmCli" "Target_NeFr40_songbirdFmFairy" "Target_NeFr40_passeri")
Runs=("Target_NeFr40" )

for x in "${Runs[@]}"; do
    echo "${x}"

    cd "${scratch_DIR}${x}_outs"

    outputDIR="${WD}/output"
    mkdir -p ${outputDIR}

    ${concat_output_PATH} elem_lik
    mv elem_lik_combined.txt ../elem_lik_combined_${x}.txt
    cp ../elem_lik_combined_${x}.txt ${outputDIR}

    ${concat_output_PATH} rate_postZ_M2
    mv rate_postZ_M2_combined.txt ../rate_postZ_M2_combined_${x}.txt
    cp ../rate_postZ_M2_combined_${x}.txt ${outputDIR}
done
# bash 3_parse_all.sh &> step3_0301.log
