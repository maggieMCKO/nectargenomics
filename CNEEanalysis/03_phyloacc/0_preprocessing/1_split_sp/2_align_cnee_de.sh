#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --ntasks=15
#SBATCH --mem=60g
#SBATCH --partition=medium
#SBATCH -o log2_alignCNEE_%A.out
#SBATCH -e log2_alignCNEE_%A.err
#SBATCH -J alignCNEE
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mcko@orn.mpg.de
##SBATCH -C "scratch|scratch2"

# adapted from https://github.com/sjswuitchik/duck_comp_gen/blob/2bf589d198122fd9ae69337a30e62ae303c7c9cb/03a_cnee_analysis/02_align_cnees.sh

module load BEDTOOLS
export cneeAna_DIR="/home/mpg08/mko/Nectar/analysis/CNEEanalysis/"

# bed
export phylop_DIR="${cneeAna_DIR}02_phylop/"
export features_bed="${phylop_DIR}galGal6_final_merged_CNEEs_named_fixchr_justchr.bed"

# maf
export maf_DIR="${cneeAna_DIR}00_maf_gff_preprocessing/2_maf/"
export input_maf="${maf_DIR}multi.anno_add1st2.maf" # added ##maf version=1 program=Bio++
export maf_processed_DIR="${maf_DIR}Processed"
# export CNEE_maf="${WD}multi.anno_CNEE.maf"

export WD="${cneeAna_DIR}03_phyloacc/preprocessing/1_way2/"
export specieslist="${WD}Sp.list"
export trans_matrix2="${cneeAna_DIR}00_maf_gff_preprocessing/galgal6a_2315.6_acc5.tsv"

export tmp_pslMappedTogalGal_dir="${WD}CneeMappedToGg6_PSL_DIR/"
export genome_fa_DIR="${cneeAna_DIR}03_phyloacc/preprocessing/0_fasta/"
export Renamed_DIR="${cneeAna_DIR}03_phyloacc/preprocessing/0_fasta/renamed_fa/"

####### RUN
### create species list [should have generated during preprocessing]
# mafSpeciesList in.maf out.lst
# mafSpeciesList ${input_maf} ${specieslist}

# 1. get ratite script from Tim, edit for specific genomes in HAL file and file names
# wget -P ${WD} https://raw.githubusercontent.com/tsackton/ratite-genomics/master/07_cnee_analysis/00_make_cnee_aligns/parse_cnee_halLiftover.pl

# 2. create a conda enviornment for running this perl script
module load conda
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
sed 's@final_cnees_long.bed@'"$features_bed"'@' ${WD}parse_cnee_halLiftover.pl | sed 's@'"$ORI_SP_set"'@'"$VAR"'@' > ${WD}parse_cnee_halLiftover_mk.pl

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
    tmpInFa="${Renamed_DIR}${SP}_renamed.fasta"
    cp ${tmpInFa} ${tmp_pslMappedTogalGal_dir}

    tmpInFaCopied="${tmp_pslMappedTogalGal_dir}${SP}_renamed.fasta"
    tmpInBed=${tmp_pslMappedTogalGal_dir}${SP}_cnees_parsed_liftover.bed
    tmpOutFa=${tmp_pslMappedTogalGal_dir}$SP.cnees.fa

    bedtools getfasta -name -s -fi ${tmpInFaCopied} -bed ${tmpInBed}  > ${tmpOutFa}

}

export -f splitsp
cut -f1 ${specieslist} | xargs -n 1 -P ${SLURM_NTASKS} -I {} bash -c 'splitsp "$@"' _ {}

# 7. use bioawk to fix up - kind of janky
rm -f "${tmp_pslMappedTogalGal_dir}all_cnees.tab"
for SP in $VAR
do

    bioawk -c fastx '{gsub(/::.*$/,"",$name); print "'"$SP"'", $name, $seq}' "${tmp_pslMappedTogalGal_dir}$SP.cnees.fa" >> "${tmp_pslMappedTogalGal_dir}all_cnees.tab"
done

# 8. check
echo "STEP8"
cut -f1,1 "${tmp_pslMappedTogalGal_dir}all_cnees.tab" | sort | uniq -c > "${tmp_pslMappedTogalGal_dir}all_cnees_summary.tab"

ls ${tmp_pslMappedTogalGal_dir} |grep "_liftover.bed" | xargs wc -l > "${tmp_pslMappedTogalGal_dir}liftover_bed_summary.tab"
