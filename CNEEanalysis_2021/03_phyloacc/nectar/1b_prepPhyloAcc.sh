#!/bin/bash
#SBATCH --time=00:05:00
#SBATCH -c 1
#SBATCH --mem=4g
#SBATCH --partition=medium
#SBATCH -o log1_NeFr40_%A.out
#SBATCH -e log1_NeFr40_%A.err
#SBATCH -J preP
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mcko@orn.mpg.de
#SBATCH -C "scratch2"

# based on
# https://github.com/sjswuitchik/duck_comp_gen/blob/master/03a_cnee_analysis/06_PhyloAcc/

# WD: /home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/03_phyloacc/oldway
source "0_setupPhyloAcc.sh"

# put named.mod, fasta, cnee bed here
Input_DIR="input/"
mkdir -p ${Input_DIR}
cp ${neutral_model_named} ${Input_DIR}
cp ${input_cnee_fa} ${Input_DIR}
cp ${input_part_bed} ${Input_DIR}

# copy to scratch
cp -r ${Input_DIR} ${scratch_DIR}

# scenario
export RUN_num="Target_NeFr40"

# make log2 dir
mkdir -p "log2_${RUN_num}/"


cd ${scratch_DIR}
# set up batch files
mkdir -p "${RUN_num}_batches" # Target_Ne40_batches
# get number of elements
num=$(wc -l ${scratch_input_part_bed} | cut -d ' ' -f 1)
# shuffle lines in random order with 0-n-1 from wc -l above
shuf -i 0-$num > ${RUN_num}_batches/full_list
cp ${RUN_num}_batches/full_list .

# make input files
# will run batches of 2000 elements each
num_per_batch=1000
split -d -a 3 -l ${num_per_batch} "${RUN_num}_batches/full_list" "${RUN_num}_batches/batch"

# set up parameter and output files
mkdir -p "${RUN_num}_param"
mkdir -p "${RUN_num}_outs"

N_array=$(( $num / ${num_per_batch} ))
N_array_m=$(( $num % ${num_per_batch} ))
if (( ${N_array_m} > 0 )) ; then
    (( N_array++ ))
    echo $N_array
fi

for I in $(seq 0 ${N_array}); # number of batches generated by line 22
do
  printf -v BATCH "%03d" $I
  PARTFILE="${scratch_DIR}${RUN_num}_batches/batch$BATCH"
  PREFIX="batch${BATCH}"
  cat > "${RUN_num}_param/run$I" <<EOF
PHYTREE_FILE ${scratch_neutral_model_named}
SEG_FILE ${scratch_input_part_bed}
ALIGN_FILE ${scratch_input_cnee_fa}
RESULT_FOLDER ${scratch_DIR}${RUN_num}_outs
PREFIX $PREFIX
ID_FILE $PARTFILE
CHAIN 1
BURNIN 1000
MCMC 7000
CONSERVE_PRIOR_A 5
CONSERVE_PRIOR_B 0.04
ACCE_PRIOR_A 10
ACCE_PRIOR_B 0.2
HYPER_GRATE_A 3
HYPER_GRATE_B 1
TARGETSPECIES HLcinPul1;HLlepAsp1;HLdicExi1;HLlicPen1;HLlicMelCas1;HLphyNov1;HLgraPic1;HLparPun1;HLlorGal1;HLtriMol2;HLaraSol1;HLamaAes1;HLphaSup1;HLfloFus1;HLcalAnn5
CONSERVE HLacaPus1;HLapuApu1;HLatrCla1;HLcalPug1;HLchaPel1;HLcliRuf1;HLcolLiv2;HLcorCor3;HLempTra1;HLfalTin1;HLfurRuf1;HLmalCya1;HLmalEle1;HLmelMel1;HLmenNov1;HLnymHol2;HLparMaj1;HLserCan1;HLstrHab1;HLtaeGut4;HLtytAlb2;aptFor1;cucCan1;falChe1;falPer1;ficAlb2;galGal6;halLeu1;opiHoa1;pseHum1
GAPCHAR -
NUM_THREAD 12
VERBOSE 0
CONSTOMIS 0.01
GAP_PROP 0.8
TRIM_GAP_PERCENT 0.8
EOF
done
