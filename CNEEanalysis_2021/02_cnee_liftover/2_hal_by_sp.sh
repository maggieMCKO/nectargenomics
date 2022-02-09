#!/bin/bash
#SBATCH --time=03:00:00 # 2h is enough for step1
##SBATCH -N 1
#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=4G # 120
#SBATCH --partition=medium
#SBATCH -o log1_hal_sp_%A.out
#SBATCH -e log1_hal_sp_%A.err
#SBATCH -J hal_sp
##SBATCH --mail-type=ALL
##SBATCH --mail-user=mcko@orn.mpg.de
#SBATCH -C "scratch2"

# WD: /home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/02_cnee_liftover

module load singularity/3.7.4 # singularity/3.2.1
export singularity_image="$HOME/Tools/cat:20200116.sif"
export specieslist="../00_inputs/sp_list.txt"
export Ref_species="galGal6"
export filtered_final="../01b_cnee_prep/bed_outputs/galGal6_final_conserved_CNEEs.bed"
export scratch_DIR="/scratch/users/$USER/phyloacc/"
mkdir -p ${scratch_DIR}

export dir_hal=liftover_proc/1_hal_sp/
export dir_psl=liftover_proc/2_psl_sp/
mkdir -p ${dir_hal} ${dir_psl}

# pigz -dc  /home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/00_inputs/multiz45way.maf.gz > ../multiz45way.maf

### RUN
# 1. comma sep species list
function join_by { local IFS="$1"; shift; echo "$*"; }
export array=($(cut -f1 ${specieslist} | grep -v "galGal6" ))
# export sp_comma=$(join_by , "${array[@]}")

# 2. maftohal
# singularity exec ${singularity_image} sh -c "cd '$HOME/Nectar/analysis/CNEEanalysis_2021/02_cnee_liftover' && maf2hal multiz45way.maf ${dir_hal}multiz45way.hal --refGenome ${Ref_species}"

# 3. halLiftover by species
runliftover () {
    export sp=$1

    singularity exec ${singularity_image} sh -c "cd '$HOME/Nectar/analysis/CNEEanalysis_2021/02_cnee_liftover' && maf2hal ../old_02_cnee_liftover/liftover_proc/0_maf/${sp}.maf ${dir_hal}${sp}.hal --refGenome ${Ref_species} && halLiftover --outPSLWithName ${dir_hal}${sp}.hal ${Ref_species} ${filtered_final} ${sp} ${dir_psl}${sp}.psl"
    # halLiftover [Options] <halFile> <srcGenome> <srcBed> <tgtGenome> <tgtBed>

}

export -f runliftover
# cut -f1 ${specieslist} | grep -v "galGal6" | parallel --max-procs ${SLURM_NTASKS} --memfree 4G --tmpdir ${scratch_DIR} runliftover {} # 3hour not enough

cut -f1 ${specieslist} | grep "HLmalEle1\|HLnymHol2\|HLparPun1\|HLstrHab1\|HLtriMol2\|HLtytAlb2\|aptFor1\|cucCan1\|falChe1\|ficAlb2\|halLeu1\|opiHoa1\|pseHum1" | parallel --max-procs ${SLURM_NTASKS} --memfree 4G --tmpdir ${scratch_DIR} runliftover {}

# 4. modify the perl script to our species set and replace the path for final cnee set
cd ${dir_psl}
export relpath_finalset="../../../01b_cnee_prep/bed_outputs/galGal6_final_conserved_CNEEs.bed"

export sp_space=$(join_by " " "${array[@]}")
echo -e "Species set: ${sp_space}"
ORI_SP_set="allMis allSin anaPla anoCar aptFor aptHaa aptOwe aptRow balReg calAnn casCas chaPel chaVoc cheMyd chrPic colLiv corBra croPor cryCin cucCan droNov eudEle falPer ficAlb fulGla gavGan halLeu lepDis melGal melUnd mesUni nipNip notPer picPub pseHum pygAde rheAme rhePen strCam taeGut tinGut"
sed 's@final_cnees_long.bed@'"$relpath_finalset"'@' ../../parse_cnee_halLiftover.pl | sed 's@'"$ORI_SP_set"'@'"$sp_space"'@' > parse_cnee_halLiftover_mk.pl

# 5. run the perl script
module load anaconda2/2019.10
#conda create -n py27 python=2.7 perl perl-app-cpanminus bedtools bioawk
source activate py27
#cpanm Math::Round # install needed perl library

perl parse_cnee_halLiftover_mk.pl # in $dir_psl


## merge
# bedtools sort -i HLcolLiv2_cnees_parsed_liftover.bed | bedtools merge -header -s -i - -d 3 -c 4 -o collapse > HLcolLiv2_cnees_parsed_liftover_merged3bp.bed # 357972 -> 357496

# bedtools sort -i aptFor1_cnees_parsed_liftover.bed | bedtools merge -header -s -i - -d 3 -c 4 -o collapse > aptFor1_cnees_parsed_liftover_merged3bp.bed

# merge () {
#     export sp=$1
#     mkdir -p ${dir_cnees}multi_mappedmerged3bp/
#     bedtools sort -i ${dir_cnees}multi_mapped/${sp}.bed  | bedtools merge -i - -d 3 > ${dir_cnees}multi_mappedmerged3bp/${sp}_merged3bp.bed
# }
# export -f merge
# cut -f1 ${specieslist} | grep -v "galGal6" | parallel --max-procs ${SLURM_NTASKS} --memfree 4G --tmpdir ${scratch_DIR} merge {}

# 6. some qc
grep "multiple_liftover_regions" final_cnees_long_liftover_parsing_log | cut -f1,1 | sort | uniq -c | sed "s/^[ \t]*//"  > multiple_liftovers_byCNEE.log

grep "no_liftover" final_cnees_long_liftover_parsing_log | cut -f1,1 | sort | uniq -c | sed "s/^[ \t]*//" > no_liftOver_byCNEE.log

find ${dir_psl} |grep "_liftover.bed" | xargs wc -l > "${dir_psl}liftover_bed_summary.tab"
