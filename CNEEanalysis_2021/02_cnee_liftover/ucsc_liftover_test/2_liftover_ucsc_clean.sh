#!/bin/bash
#SBATCH --time=03:30:00 # 1:30
#SBATCH -N 1
#SBATCH --ntasks=10
#SBATCH --mem=120G # 120
#SBATCH --partition=medium
#SBATCH -o log1_liftover_%A.out
#SBATCH -e log1_liftover_%A.err
#SBATCH -J liftover
##SBATCH --mail-type=ALL
##SBATCH --mail-user=mcko@orn.mpg.de
#SBATCH -C "scratch2"

source "0_setup.sh"

mkdir -p Sp_psl Sp_chain Sp_chainMerge Sp_sort Sp_2bit Sp_chromSizes Sp_net Sp_liftover

### RUN
# cd ${maf_DIR}
${ucsc_path}faToTwoBit ../00_inputs/genomes/galGal6.fa Sp_2bit/galGal6.2bit
${ucsc_path}twoBitInfo Sp_2bit/galGal6.2bit stdout | sort -k2,2nr > Sp_chromSizes/galGal6.chrom.sizes

runliftover () {
    export sp=$1

    ## get and .2bit chrom.sizes
    ${ucsc_path}faToTwoBit ${tNibDir} Sp_2bit/$sp.2bit
    ${ucsc_path}twoBitInfo Sp_2bit/$sp.2bit stdout | sort -k2,2nr > Sp_chromSizes/$sp.chrom.sizes

    ## maf to psl
    maf_input="Sp_maf/${sp}.maf"
    psl_output="Sp_psl/${sp}.psl"
    "${last_path}maf-convert" psl ${maf_input} > ${psl_output}

    ## Chaining - join clost alignments
    tNibDir=$(find ../00_inputs/genomes/ | grep ${sp})
    qNibDir="../00_inputs/genomes/galGal6.fa"
    # axtChain [options] -linearGap=loose in.axt tNibDir qNibDir out.chain
    ${ucsc_path}axtChain -linearGap=medium -psl -faT -faQ ${psl_output} ${tNibDir} ${qNibDir} Sp_chain/t${sp}_qgalGal6.chain

    # tNibDir=$(find Sp_2bit | grep ${sp})
    # qNibDir=$(find Sp_2bit | grep "galGal6")
    # ${ucsc_path}axtChain -linearGap=medium -psl ${psl_output} ${tNibDir} ${qNibDir} Sp_chain/t${sp}_qgalGal6.chain

    # merge and split
    MYSCRATCH=`mktemp -d ${scratch_DIR}tmp.XXXXXXXX`
    mkdir -p Sp_chainMerge/$sp
    ${ucsc_path}chainMergeSort -saveId -tempDir=${MYSCRATCH} Sp_chain/t${sp}_qgalGal6.chain | ${ucsc_path}chainSplit Sp_chainMerge/$sp stdin -lump=50
    rm -rf ${MYSCRATCH}

    # concat and sort
    mkdir -p Sp_sort/$sp
    cat Sp_chainMerge/$sp/*.chain > Sp_sort/$sp/$sp.chain
    ${ucsc_path}chainSort Sp_sort/$sp/$sp.chain Sp_sort/$sp/$sp_sorted.chain

    ## Netting
    # chainNet in.chain target.sizes query.sizes target.net query.net
    ${ucsc_path}chainNet Sp_sort/$sp/$sp_sorted.chain Sp_chromSizes/$sp.chrom.sizes Sp_chromSizes/galGal6.chrom.sizes Sp_net/$sp.net /dev/null

    # netChainSubset in.net in.chain out.chain
    ${ucsc_path}netChainSubset Sp_net/$sp.net Sp_sort/$sp/$sp_sorted.chain Sp_liftover/galGal6_to_${sp}.liftOver

}

export -f runliftover

cut -f1 ${specieslist} | grep -v "galGal6" | parallel --max-procs 4 --memfree 4G --tmpdir ${scratch_DIR} runliftover {}
