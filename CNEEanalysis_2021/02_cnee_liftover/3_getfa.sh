#!/bin/bash
#SBATCH --time=01:00:00 # 2h is enough for step1
##SBATCH -N 1
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=5G # 120
#SBATCH --partition=medium
#SBATCH -o log3_getfa_%A.out
#SBATCH -e log3_getfa_%A.err
#SBATCH -J getfa
##SBATCH --mail-type=ALL
##SBATCH --mail-user=mcko@orn.mpg.de
#SBATCH -C "scratch2"

# WD: /home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/02_cnee_liftover
source "0_setup.sh"

export dir_hal=liftover_proc/1_hal_sp/
export dir_psl=liftover_proc/2_psl_sp/
export dir_cneefa=liftover_proc/3_fa_sp/
export dir_cneeconcatfa=liftover_proc/4_concatfa/
mkdir -p ${dir_hal} ${dir_psl} ${dir_cneefa} ${dir_cneeconcatfa}

### RUN

getfa (){
    export sp=$1

    qfa=$(find ../00_inputs/genomes/ | grep "fa$"| grep ${sp})
    tmpOutFa=${dir_cneefa}${sp}_cnees.fa

    if [ "$sp" == "galGal6" ]; then
        bedtools getfasta -name -s -fi ${qfa} -bed ${filtered_final}  > ${tmpOutFa}
    else
        bedtools getfasta -name -s -fi ${qfa} -bed ${dir_psl}${sp}_cnees_parsed_liftover.bed  > ${tmpOutFa}
    fi
}
export -f getfa
cut -f1 ${specieslist} | parallel --max-procs ${SLURM_NTASKS} --memfree 4G --tmpdir ${scratch_DIR} getfa {}
# ~20 mins

# adapted from https://github.com/sjswuitchik/duck_comp_gen/blob/2bf589d198122fd9ae69337a30e62ae303c7c9cb/03a_cnee_analysis/02_align_cnees.sh

## use bioawk to fix up - kind of janky
# 1 prepare the species set
function join_by { local IFS="$1"; shift; echo "$*"; }
export array=($(cut -f1 ${specieslist} | grep -v "galGal6" ))
export sp_space=$(join_by " " "${array[@]}")
echo -e "Species set: ${sp_space}"

module load anaconda2/2019.10
#conda create -n py27 python=2.7 perl perl-app-cpanminus bedtools bioawk
source activate py27
#cpanm Math::Round # install needed perl library


rm -f "${dir_cneeconcatfa}all_cnees.tab"
for SP in $sp_space
do
    bioawk -c fastx '{gsub(/::.*$/,"",$name); print "'"$SP"'", $name, $seq}'  "${dir_cneefa}${SP}_cnees.fa" >> "${dir_cneeconcatfa}all_cnees.tab"
done

cut -f1,1 "${dir_cneeconcatfa}all_cnees.tab" | sort | uniq -c > "${dir_cneeconcatfa}all_cnees_summary.tab"
