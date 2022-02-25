#!/bin/bash
#SBATCH --time=02:30:00 # 2h is enough for step1
##SBATCH -N 1
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=5G # 120
#SBATCH --partition=medium
#SBATCH -o log1_liftover_%A.out
#SBATCH -e log1_liftover_%A.err
#SBATCH -J liftover
##SBATCH --mail-type=ALL
##SBATCH --mail-user=mcko@orn.mpg.de
#SBATCH -C "scratch2"

source "0_setup.sh"

export dir_2bit=liftover_proc/1_2bit/
export dir_chromsizes=liftover_proc/2_chromSizes/
export dir_psl=liftover_proc/3_psl/
export dir_chain=liftover_proc/4_chain/
export dir_chainMerged=liftover_proc/5_chainMerge/
export dir_chainsorted=liftover_proc/6_chainsorted/
export dir_chainMerged2=liftover_proc/7_chainsortedMerged/
export dir_net=liftover_proc/8_net/
export dir_liftover=liftover_proc/9_liftover/
export dir_cnees=liftover_proc/10_cnees/
export dir_cneefa=liftover_proc/11_fa/
export dir_cneeconcatfa=liftover_proc/12_concatfa/
mkdir -p liftover_proc ${dir_2bit} ${dir_chromsizes} ${dir_psl} ${dir_chain} ${dir_chainMerged} ${dir_chainsorted} ${dir_chainMerged2} ${dir_net} ${dir_liftover} ${dir_cnees} ${dir_cneefa} ${dir_cneeconcatfa}

### RUN
# cd ${maf_DIR}
${ucsc_path}faToTwoBit ../00_inputs/genomes/galGal6.fa ${dir_2bit}galGal6.2bit
${ucsc_path}twoBitInfo ${dir_2bit}galGal6.2bit stdout | sort -k2,2nr > ${dir_chromsizes}galGal6.chrom.sizes

runliftover () {
    export sp=$1

    # for "liftOver", the map.chain file has the old genome as the target and the new genome as the query.

    ## get chrom.sizes
    qfa=$(find ../00_inputs/genomes/ | grep "fa$"| grep ${sp})
    ${ucsc_path}faToTwoBit ${qfa} ${dir_2bit}${sp}.2bit
    ${ucsc_path}twoBitInfo ${dir_2bit}${sp}.2bit stdout | sort -k2,2nr > ${dir_chromsizes}${sp}.chrom.sizes

    ## maf to psl
    maf_input="liftover_proc/0_maf/${sp}.maf"
    psl_output="${dir_psl}${sp}.psl"
    mafToPsl ${sp} ${Ref_species} ${maf_input} ${psl_output}
    # mafToPsl querySrc targetSrc in.maf out.psl

    ## Chaining - join clost alignments
    tNibDir="${dir_2bit}galGal6.2bit"
    qNibDir=$(find ${dir_2bit} | grep ${sp})
    ${ucsc_path}axtChain -linearGap=medium -psl ${psl_output} ${tNibDir} ${qNibDir} ${dir_chain}tgalGal6_q${sp}.chain
    # axtChain [options] -linearGap=loose in.axt tNibDir qNibDir out.chain

    ## merge and split
    MYSCRATCH=`mktemp -d ${scratch_DIR}tmp.XXXXXXXX`
    mkdir -p ${dir_chainMerged}${sp}
    # chainMergeSort file(s)
    # chainSplit outDir inChain(s)
    ${ucsc_path}chainMergeSort -saveId -tempDir=${MYSCRATCH} ${dir_chain}tgalGal6_q${sp}.chain | ${ucsc_path}chainSplit ${dir_chainMerged}${sp} stdin -lump=50

    ## sort
    mkdir -p ${dir_chainsorted}${sp}
    for i in ${dir_chainMerged}${sp}/*.chain; do ${ucsc_path}chainSort $i ${dir_chainsorted}${sp}/`basename $i`; done
    # chainSort inFile outFile

    ${ucsc_path}chainMergeSort -saveId -tempDir=${MYSCRATCH} ${dir_chainsorted}${sp}/*.chain > ${dir_chainMerged2}${sp}_MergedSorted.chain
    rm -rf ${MYSCRATCH}
    # saveId or not seems the same

    ## Netting
    # chainNet in.chain target.sizes query.sizes target.net query.net
    ${ucsc_path}chainNet ${dir_chainMerged2}${sp}_MergedSorted.chain ${dir_chromsizes}galGal6.chrom.sizes ${dir_chromsizes}${sp}.chrom.sizes ${dir_net}t_${sp}.net ${dir_net}q_${sp}.net  # /dev/null

    # netChainSubset in.net in.chain out.chain
    ${ucsc_path}netChainSubset ${dir_net}t_${sp}.net ${dir_chainMerged2}${sp}_MergedSorted.chain ${dir_liftover}galGal6_to_${sp}.liftOver

    # liftover
    # liftOver oldFile map.chain newFile unMapped
    mkdir -p ${dir_cnees}mapped ${dir_cnees}unmapped
    ${ucsc_path}liftOver ${filtered_final} ${dir_liftover}galGal6_to_${sp}.liftOver ${dir_cnees}mapped/${sp}.bed ${dir_cnees}unmapped/${sp}.bed
    # The map.chain file has the old genome as the target and the new genome as the query.
}

export -f runliftover
cut -f1 ${specieslist} | grep -v "galGal6" | parallel --max-procs ${SLURM_NTASKS} --memfree 4G --tmpdir ${scratch_DIR} runliftover {}


merge () {
    export sp=$1
    mkdir -p ${dir_cnees}mappedmerged3bp/
    bedtools sort -i ${dir_cnees}mapped/${sp}.bed  | bedtools merge -i - -d 3 > ${dir_cnees}mappedmerged3bp/${sp}_merged3bp.bed
}
export -f merge
cut -f1 ${specieslist} | grep -v "galGal6" | parallel --max-procs ${SLURM_NTASKS} --memfree 4G --tmpdir ${scratch_DIR} merge {}


getfa (){
    export sp=$1

    qfa=$(find ../00_inputs/genomes/ | grep "fa$"| grep ${sp})
    tmpOutFa=${dir_cneefa}${sp}_cnees.fa

    if [ "$sp" == "galGal6" ]; then
        bedtools getfasta -name -s -fi ${qfa} -bed ${filtered_final}  > ${tmpOutFa}
    else
        bedtools getfasta -name -s -fi ${qfa} -bed ${dir_cnees}mapped/${sp}.bed  > ${tmpOutFa}
    fi
}
export -f getfa
cut -f1 ${specieslist} | parallel --max-procs ${SLURM_NTASKS} --memfree 4G --tmpdir ${scratch_DIR} getfa {}
# ~20 mins

# adapted from https://github.com/sjswuitchik/duck_comp_gen/blob/2bf589d198122fd9ae69337a30e62ae303c7c9cb/03a_cnee_analysis/02_align_cnees.sh

## use bioawk to fix up - kind of janky
# rm -f "${dir_cneeconcatfa}all_cnees.tab"
# 1 prepare the species set
SP_set=$(cut -f1 ${specieslist})
VAR=""
for ELEMENT in ${SP_set}; do
  VAR+="${ELEMENT} "
done
echo -e "Species set: ${VAR}"

module load anaconda2/2019.10
#conda create -n py27 python=2.7 perl perl-app-cpanminus bedtools bioawk
source activate py27
#cpanm Math::Round # install needed perl library

for SP in $VAR
do
    bioawk -c fastx '{gsub(/::.*$/,"",$name); print "'"$SP"'", $name, $seq}'  "${dir_cneefa}${SP}_cnees.fa" >> "${dir_cneeconcatfa}all_cnees.tab"
done

cut -f1,1 "${dir_cneeconcatfa}all_cnees.tab" | sort | uniq -c > "${dir_cneeconcatfa}all_cnees_summary.tab"
