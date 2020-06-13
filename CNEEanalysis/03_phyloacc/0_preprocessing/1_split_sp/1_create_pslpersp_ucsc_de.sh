#!/bin/bash
#SBATCH --time=01:30:00
#SBATCH --ntasks=30
#SBATCH --mem=120g
#SBATCH --partition=medium
#SBATCH -o log1_maf2psl_%A.out
#SBATCH -e log1_maf2psl_%A.err
#SBATCH -J maf2psl
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mcko@orn.mpg.de
#SBATCH -C "scratch|scratch2"

module load UCSC
module load BEDTOOLS

# ls /usr/product/bioinfo/SL_7.0/BIOINFORMATICS/UCSC/20160601 | grep -i psl
export cneeAna_DIR="/home/mpg08/mko/Nectar/analysis/CNEEanalysis/"

# bed
export phylop_DIR="${cneeAna_DIR}02_phylop/"
export features_bed="${phylop_DIR}galGal6_final_merged_CNEEs_named_fixchr_justchr.bed"

# maf
export maf_DIR="${cneeAna_DIR}00_maf_gff_preprocessing/2_maf/"
export input_maf="${maf_DIR}multi.anno_add1st2.maf" # added ##maf version=1 program=Bio++
export maf_processed_DIR="${maf_DIR}Processed"
# export CNEE_maf="${WD}multi.anno_CNEE.maf"

export WD="${cneeAna_DIR}03_phyloacc/preprocessing/1_split_sp/"
export specieslist="${WD}Sp.list"
export trans_matrix2="${cneeAna_DIR}00_maf_gff_preprocessing/galgal6a_2315.6_acc5.tsv"

export cneePSL="${WD}cnee.psl"
### create species list
# mafSpeciesList in.maf out.lst
mafSpeciesList ${input_maf} ${specieslist}

### convert cnee.bed to cnee.psl
awk '$3 ~ /region/ {print $1 "\t" $5}' ${cneeAna_DIR}GCF_000002315.6_GRCg6a_genomic_chrtabFixed.gff > ${WD}galGal6chromSize.txt
bedToPsl ${WD}galGal6chromSize.txt ${features_bed} ${cneePSL} # works



splitChr (){
    Chr=$1
    tmpSP=$2
    echo -e "tmpSP: ${tmpSP}"

    ## input maf ##
    chr_maf="${Chr}.maf"
    tmp_maf=$(find ${maf_processed_DIR} | grep ${chr_maf})
    echo -e "tmp_maf: ${tmp_maf}"
    line_maf=$(wc -l $tmp_maf)
    echo -e "maf: ${line_maf}"
    line_maf_cut=$(echo ${line_maf} | cut -d ' ' -f1 )

    if (( $line_maf_cut > 2 )) ; then
        echo $Chr

        ## convert species maf to psl
        tmp_sp_psl_dir="${WD}tmp_PSL_DIR/${tmpSP}_PSL_DIR/"
        mkdir -p ${tmp_sp_psl_dir}
        tmp_psl="${tmp_sp_psl_dir}${tmpSP}_${Chr}.psl"
        echo -e "tmp_psl: ${tmp_psl}"
        mafToPsl ${tmpSP} galGal6  ${tmp_maf} ${tmp_psl}

    fi
}

export -f splitChr




splitsp (){
    export SP=$1

    ### exec splitchr
    cut -f2 ${trans_matrix2} | xargs -n 1 -P 10 -I {} bash -c 'splitChr "$@" ${SP}' _ {}

    ### concate/sort chr.psl
    tmp_sp_psl_dir="${WD}tmp_PSL_DIR/${SP}_PSL_DIR/" # from
    tmp_psl_dir="${WD}ConCatSort_PSL_DIR/"                      # to
    mkdir -p ${tmp_psl_dir}
    tmp_psl="${tmp_psl_dir}${SP}.psl" # important: in species.psl format for Tim's perl script

    echo -e "tmp_psl: ${tmp_psl}"
    echo -e "MYSCRATCH: ${MYSCRATCH}"
    echo -e "tmp_sp_psl_dir: ${tmp_sp_psl_dir}"
    mkdir -p /scratch/${USER}
    MYSCRATCH=`mktemp -d /scratch/${USER}/tmp.XXXXXXXX`
    pslSort -nohead dirs ${tmp_psl} ${MYSCRATCH} ${tmp_sp_psl_dir} # [ConCatSort_PSL_DIR]
    # pslSort dirs[1|2] outFile tempDir inDir(s)OrFile(s)
    rm -rf ${MYSCRATCH}

    # map cnee
    tmp_pslMapped_dir="${WD}CneeMapped_PSL_DIR/"                      # to
    mkdir -p ${tmp_pslMapped_dir}
    tmp_pslMapped="${tmp_pslMapped_dir}${SP}.psl" # important: in species.psl format for Tim's perl script

    tmp_pslMapInfo_dir="${WD}CneeMapInfo_PSL_DIR/"                     # to
    mkdir -p ${tmp_pslMapInfo_dir}
    tmp_pslMapInfo="${tmp_pslMapInfo_dir}${SP}_mapInfo"
    pslMap -swapMap -mapInfo=${tmp_pslMapInfo} ${tmp_psl} ${cneePSL} ${tmp_pslMapped}

    # add CneeID to queryname
    tmp_tmp_pslMappedmod_dir="${WD}CneeMappedMod_PSL_DIR/"
    mkdir -p ${tmp_tmp_pslMappedmod_dir}
    tmp_mod1="${tmp_tmp_pslMappedmod_dir}${SP}_1.psl"
    awk -F "\t" '{OFS=FS}{$10=$10"__"$14; print }' ${tmp_pslMapped} > ${tmp_mod1} # double underscore

    # re-map to galGal6 corrodinates
    tmp_pslMappedToGg6="${tmp_tmp_pslMappedmod_dir}${SP}_2_remappedToGgal6.psl"
    tmp_pslMapInfoToGg6="${tmp_tmp_pslMappedmod_dir}${SP}_2_remappedToGgal6_mapInfo"
    pslMap -mapInfo=${tmp_pslMapInfoToGg6} ${tmp_mod1} ${cneePSL} ${tmp_pslMapInfoToGg6}

    # swap query and target
    tmp_mod3="${tmp_tmp_pslMappedmod_dir}${SP}_3.psl"
    pslSwap ${tmp_pslMapInfoToGg6} ${tmp_mod3}

    # sort by query
    mkdir -p /scratch/${USER}
    MYSCRATCH=`mktemp -d /scratch/${USER}/tmp.XXXXXXXX`
    tmp_mod4="${tmp_tmp_pslMappedmod_dir}${SP}_4_sort.psl"
    pslSort -nohead dirs ${tmp_mod4} ${MYSCRATCH} ${tmp_mod3}
    # pslSort dirs[1|2] outFile tempDir inDir(s)OrFile(s)
    rm -rf ${MYSCRATCH}

    # break up target name (targetChr_CneeID), and move CneeID to the 1st col, and fix strand
    tmp_pslMappedTogalGal_dir="${WD}CneeMappedToGg6_PSL_DIR/"
    mkdir -p ${tmp_pslMappedTogalGal_dir}
    tmp_mod5="${tmp_pslMappedTogalGal_dir}${SP}.psl" # important: in species.psl format for Tim's perl script
    awk '{split($0,a,"__"); print a[1], a[2]}' OFS="\t" ${tmp_mod4} | awk '{ print $15,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$16,$17,$18,$19,$20,$21,$22}' OFS="\t" | awk -F "\t" '{OFS=FS}{ $10="+"$10 ; print  }' > ${tmp_mod5} # double underscore

}

export -f splitsp

# cut -f1 ${specieslist} | xargs -n 1 -P 10 -I {} bash -c 'splitsp "$@"' _ {}
cut -f1 ${specieslist} | xargs -n 1 -P ${SLURM_NTASKS} -I {} bash -c 'splitsp "$@"' _ {}
