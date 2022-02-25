#!/bin/bash
#SBATCH --time=03:30:00 # 1:30
#SBATCH -N 1
#SBATCH --ntasks=24
#SBATCH --mem=400 # 120
#SBATCH --partition=fat
#SBATCH -o log1_maf2psl_%A.out
#SBATCH -e log1_maf2psl_%A.err
#SBATCH -J maf2psl
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mcko@orn.mpg.de
#SBATCH -C "scratch2"

source "0_setup.sh"

# use scratch: setup in 0_setup.sh

### convert cnee.bed to cnee.psl [run separately works]
# awk '$3 ~ /region/ {print $1 "\t" $5}' ${fixed_chicken_anno} | awk '{gsub(/galGal6./, "") ; print}' > ${WD}galGal6chromSize.txt # (.chr) # works
# bedToPsl ${WD}galGal6chromSize.txt ${filtered_final} ${cneePSL} # works

### to extract per Chr. maf of species X to psl on galgal6 cord.
mafToPslbyChr (){
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
        tmp_sp_psl_dir="${WD}d1_tmp_PSL_DIR/${tmpSP}_PSL_DIR/"
        mkdir -p ${tmp_sp_psl_dir}
        tmp_psl="${tmp_sp_psl_dir}${tmpSP}_${Chr}.psl"
        echo -e "tmp_psl: ${tmp_psl}"
        mafToPsl ${tmpSP} ${Ref_species} ${tmp_maf} ${tmp_psl}
        # mafToPsl querySrc targetSrc in.maf out.psl
    fi
}
export -f mafToPslbyChr


splitsp (){
    export SP=$1

    ### create all subfolders
    tmp_psl_dir="${WD}d2_ConCatSort_PSL_DIR/"           # to
    mkdir -p ${tmp_psl_dir}
    tmp_pslMapped_dir="${WD}d3_1_CneeMapped_PSL_DIR/"   # to
    mkdir -p ${tmp_pslMapped_dir}
    # tmp_pslMapInfo_dir="${WD}d3_2_CneeMapInfo_PSL_DIR/"  # to
    # mkdir -p ${tmp_pslMapInfo_dir}
    tmp_tmp_pslMappedmod_dir="${WD}d4_CneeMappedMod_PSL_DIR/"
    mkdir -p ${tmp_tmp_pslMappedmod_dir}
    tmp_pslMappedTogalGal_dir="${WD}d5_CneeMappedToGg6_PSL_DIR/"
    mkdir -p ${tmp_pslMappedTogalGal_dir}
    # final output
    tmp_mod5="${tmp_pslMappedTogalGal_dir}${SP}.psl"
    # important: in species.psl format (q:galgal; t: sp) for Tim's perl script

    if [[ -s ${tmp_mod5} ]]
	then
		# file exists and is not empty
		:
	else

        ### exec mafToPslbyChr: to extract per Chr. maf of species X to psl on galgal6 cord.
        cut -f2 ${trans_matrix_for_maf} | xargs -n 1 -P 10 -I {} bash -c 'mafToPslbyChr "$@" ${SP}' _ {}

        ### concate/sort chr.psl
        tmp_sp_psl_dir="${WD}d1_tmp_PSL_DIR/${SP}_PSL_DIR/" # from
        tmp_psl="${tmp_psl_dir}${SP}.psl"
        # important: in species.psl format for Tim's perl script
        echo -e "tmp_psl: ${tmp_psl}"

        MYSCRATCH=`mktemp -d ${scratch_DIR}tmp.XXXXXXXX`
        echo -e "MYSCRATCH: ${MYSCRATCH}"
        pslSort -nohead dirs ${tmp_psl} ${MYSCRATCH} ${tmp_sp_psl_dir} # [ConCatSort_PSL_DIR]
        # pslSort: Merge and sort psCluster .psl output files
        # pslSort dirs[1|2] outFile tempDir inDir(s)OrFile(s)
        rm -rf ${MYSCRATCH}

        ### map cnees to species
        tmp_pslMapped="${tmp_pslMapped_dir}${SP}.psl"
        tmp_pslMapInfo="${tmp_pslMapInfo_dir}${SP}_mapInfo"
        # pslMap -swapMap -mapInfo=${tmp_pslMapInfo} ${tmp_psl} ${cneePSL} ${tmp_pslMapped}
        pslMap -swapMap ${tmp_psl} ${cneePSL} ${tmp_pslMapped}
        # pslMap: map PSLs alignments to new targets
        # pslMap [options] inPsl mapFile outPsl
        # Given inPsl and mapPsl, where the target of inPsl is the query of mapPsl
        # tmp_psl: q=sp; t=galgal
        # cneePSL: q=cnee; t=galgal
        # -swapMap - swap query and target sides of map file.
        # tmp_pslMapped: q=sp; t=cnee (cnees mapped to species)

        # add CneeID to queryname
        tmp_mod1="${tmp_tmp_pslMappedmod_dir}${SP}_1.psl"
        awk -F "\t" '{OFS=FS}{$10=$10"__"$14; print }' ${tmp_pslMapped} > ${tmp_mod1} # double underscore

        # re-map to galGal6 corrodinates
        tmp_pslMappedToGg6="${tmp_tmp_pslMappedmod_dir}${SP}_2_remappedToGgal6.psl"
        # tmp_pslMapInfoToGg6="${tmp_tmp_pslMappedmod_dir}${SP}_2_remappedToGgal6_mapInfo"
        # pslMap -mapInfo=${tmp_pslMapInfoToGg6} ${tmp_mod1} ${cneePSL} ${tmp_pslMapInfoToGg6}
        pslMap ${tmp_mod1} ${cneePSL} ${tmp_pslMapInfoToGg6}
        # tmp_mod1: q=sp; t=cnee
        # cneePSL: q=cnee; t=galgal
        # tmp_pslMapInfoToGg6: q=sp; t=galgal (cnees mapped to species and liftover to ggal6)

        # swap query and target
        tmp_mod3="${tmp_tmp_pslMappedmod_dir}${SP}_3.psl"
        pslSwap ${tmp_pslMapInfoToGg6} ${tmp_mod3}
        # tmp_mod3: q=galgal; t=sp (cnee)

        # sort by query
        MYSCRATCH=`mktemp -d ${scratch_DIR}tmp.XXXXXXXX`
        echo -e "MYSCRATCH: ${MYSCRATCH}"
        tmp_mod4="${tmp_tmp_pslMappedmod_dir}${SP}_4_sort.psl"
        pslSort -nohead dirs ${tmp_mod4} ${MYSCRATCH} ${tmp_mod3}
        # pslSort outFile tempDir inDir(s)OrFile(s)
        # tmp_mod4: q=galgal; t=sp (cnee sorted by gaggal6)
        rm -rf ${MYSCRATCH}

        # break up target name (targetChr_CneeID), and move CneeID to the 1st col, and fix strand
        awk '{split($0,a,"__"); print a[1], a[2]}' OFS="\t" ${tmp_mod4} | awk '{ print $15,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$16,$17,$18,$19,$20,$21,$22}' OFS="\t" | awk -F "\t" '{OFS=FS}{ $10="+"$10 ; print  }' > ${tmp_mod5} # double underscore

	fi
}
export -f splitsp

cut -f1 ${specieslist} | parallel --memfree 14 --tmpdir ${scratch_DIR} splitsp {}
