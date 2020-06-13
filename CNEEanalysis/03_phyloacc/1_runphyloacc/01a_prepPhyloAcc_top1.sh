#!/bin/bash
#SBATCH --time=00:05:00
#SBATCH -c 1
#SBATCH --mem=1g
#SBATCH --partition=medium
#SBATCH -o log1_preP_%A.out
#SBATCH -e log1_preP_%A.err
#SBATCH -J preP
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mcko@orn.mpg.de
#SBATCH -C "scratch2"

# based on
# https://github.com/sjswuitchik/duck_comp_gen/blob/master/03a_cnee_analysis/06_PhyloAcc/

source /home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/1_run_phyloacc/setupPhyloAcc.sh

# put named.mod, fasta, cnee bed here
Input_DIR="${WD}input"
mkdir -p ${Input_DIR}
cp ${neutral_model_named} ${Input_DIR}
cp ${input_cnee_fa} ${Input_DIR}
cp ${input_part_bed} ${Input_DIR}

# copy to scratch
cp -r ${Input_DIR} ${MYSCRATCH}
cd ${MYSCRATCH}


RUN_num="top1"

# set up batch files
mkdir -p "${RUN_num}_batches" # top1_batches
# get number of elements
num=$(wc -l ${scartch_input_part_bed} | cut -d ' ' -f 1)
# shuffle lines in random order with 0-n-1 from wc -l above
shuf -i 0-$num > ${RUN_num}_batches/full_list
cp ${RUN_num}_batches/full_list ${WD}

# make input files
# will run batches of 2000 elements each
num_per_batch=2000
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
  PARTFILE="${MYSCRATCH}${RUN_num}_batches/batch$BATCH"
  PREFIX="batch${BATCH}"
  cat > "${RUN_num}_param/run$I" <<EOF
PHYTREE_FILE ${scartch_neutral_model_named}
SEG_FILE ${scartch_input_part_bed}
ALIGN_FILE ${scartch_input_cnee_fa}
RESULT_FOLDER ${MYSCRATCH}${RUN_num}_outs
PREFIX $PREFIX
ID_FILE $PARTFILE
CHAIN 1
BURNIN 1000
MCMC 4000
CONSERVE_PRIOR_A 5
CONSERVE_PRIOR_B 0.04
ACCE_PRIOR_A 10
ACCE_PRIOR_B 0.2
HYPER_GRATE_A 3
HYPER_GRATE_B 1
TARGETSPECIES HLphyNov1;HLlicCas1;HLtriMol2;HLcalAnn5;HLfloFus1
CONSERVE HLacaPus1;HLcalPug1;HLchaPel1;HLcolLiv2;HLcorCor3;HLempTra1;HLfalTin1;HLfurRuf1;HLmalCya1;HLparMaj1;HLserCan1;HLstrHab1;HLtaeGut4;HLtytAlb2;HLzosLat1;aptFor1;cucCan1;falChe1;falPer1;ficAlb2;galGal6;halLeu1;melUnd1;opiHoa1;pseHum1
GAPCHAR -
NUM_THREAD 12
VERBOSE 0
CONSTOMIS 0.01
GAP_PROP 0.8
TRIM_GAP_PERCENT 0.8
EOF
done
