#!/bin/bash
#SBATCH --time=01:00:00 # only use 26 mins 20210819
#SBATCH --ntasks=15
#SBATCH --mem=60g
#SBATCH --partition=medium
#SBATCH -o log2_alignCNEE_%A.out
#SBATCH -e log2_alignCNEE_%A.err
#SBATCH -J alignCNEE
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mcko@orn.mpg.de
#SBATCH -C "scratch|scratch2"

# adapted from https://github.com/sjswuitchik/duck_comp_gen/blob/2bf589d198122fd9ae69337a30e62ae303c7c9cb/03a_cnee_analysis/02_align_cnees.sh

source "0_setup.sh"



####### RUN
# 1. get ratite script from Tim, edit for specific genomes in HAL file and file names
# wget -P ${WD} https://raw.githubusercontent.com/tsackton/ratite-genomics/master/07_cnee_analysis/00_make_cnee_aligns/parse_cnee_halLiftover.pl

# 2. create a conda enviornment for running this perl script
module load anaconda2/2019.10
#conda create -n py27 python=2.7 perl perl-app-cpanminus bedtools bioawk
source activate py27
#cpanm Math::Round # install needed perl library

# 3. modify the perl script to my purpose
# 3.1 prepare the species set
SP_set=$(cut -f1 ${specieslist})
VAR=""
for ELEMENT in ${SP_set}; do
  VAR+="${ELEMENT} "
done
echo -e "Species set: ${VAR}"

# 3.2 modify
ORI_SP_set="allMis allSin anaPla anoCar aptFor aptHaa aptOwe aptRow balReg calAnn casCas chaPel chaVoc cheMyd chrPic colLiv corBra croPor cryCin cucCan droNov eudEle falPer ficAlb fulGla gavGan halLeu lepDis melGal melUnd mesUni nipNip notPer picPub pseHum pygAde rheAme rhePen strCam taeGut tinGut"
sed 's@final_cnees_long.bed@'"$filtered_final"'@' ${WD}parse_cnee_halLiftover.pl | sed 's@'"$ORI_SP_set"'@'"$VAR"'@' > ${WD}parse_cnee_halLiftover_mk.pl

# # 4. run the perl script
cd ${tmp_pslMappedTogalGal_dir}
perl ${WD}parse_cnee_halLiftover_mk.pl

# 5. some qc
grep "multiple_liftover_regions" ${tmp_pslMappedTogalGal_dir}final_cnees_long_liftover_parsing_log | cut -f1,1 | sort | uniq -c | sed "s/^[ \t]*//"  > ${tmp_pslMappedTogalGal_dir}multiple_liftovers_byCNEE.log

grep "no_liftover" ${tmp_pslMappedTogalGal_dir}final_cnees_long_liftover_parsing_log | cut -f1,1 | sort | uniq -c | sed "s/^[ \t]*//" > ${tmp_pslMappedTogalGal_dir}no_liftOver_byCNEE.log

# 6. generating mafft input: bedtools to get a fasta file for each species with an entry for each CNEE, then split and merge by CNEE name, output should be one file per CNEE with fasta header = species (instead of one file per species with fasta header = CNEE)

splitsp (){
    SP=$1

    # # for some reason, bedtools getfasta needs to be exe-ed in the dir of in_fa, in_bed, and out_fa
    # tmpInFa="${Renamed_DIR}${SP}_renamed.fasta"
    tmpInFa="${Renamed_DIR}${SP}.fa" # 2021
    cp ${tmpInFa} ${tmp_pslMappedTogalGal_dir}

    tmpInFaCopied="${tmp_pslMappedTogalGal_dir}${SP}.fa"
    tmpInBed=${tmp_pslMappedTogalGal_dir}${SP}_cnees_parsed_liftover.bed
    tmpOutFa=${tmp_pslMappedTogalGal_dir}$SP.cnees.fa

    bedtools getfasta -name -s -fi ${tmpInFaCopied} -bed ${tmpInBed}  > ${tmpOutFa}

}

export -f splitsp
# cut -f1 ${specieslist} | xargs -n 1 -P ${SLURM_NTASKS} -I {} bash -c 'splitsp "$@"' _ {}
cut -f1 ${specieslist} | parallel splitsp {}

# 7. use bioawk to fix up - kind of janky
# rm -f "${tmp_pslMappedTogalGal_dir}all_cnees.tab"
for SP in $VAR
do

    bioawk -c fastx '{gsub(/::.*$/,"",$name); print "'"$SP"'", $name, $seq}' "${tmp_pslMappedTogalGal_dir}$SP.cnees.fa" >> "${tmp_pslMappedTogalGal_dir}all_cnees.tab"
done

# 8. check
echo "STEP8"
cut -f1,1 "${tmp_pslMappedTogalGal_dir}all_cnees.tab" | sort | uniq -c > "${tmp_pslMappedTogalGal_dir}all_cnees_summary.tab"

ls ${tmp_pslMappedTogalGal_dir} |grep "_liftover.bed" | xargs wc -l > "${tmp_pslMappedTogalGal_dir}liftover_bed_summary.tab"
