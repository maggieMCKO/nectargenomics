#!/bin/bash

source setup_paths.sh

runPhyloPbyChr (){
    Chr=$1
    chr_clean=$(echo $Chr | sed 's/galGal6.//')
    # echo $Chr

    # input maf
    chr_maf="${chr_clean}.maf"
    tmp_maf=$(find ${maf_processed_DIR} | grep ${chr_maf})
    line_maf=$(wc -l $tmp_maf)
    echo -e "maf: ${line_maf}"
    line_maf_cut=$(echo ${line_maf} | cut -d ' ' -f1 )

    # input bed
    tmp_bed="${CNEEbyChr}galGal6_${chr_clean}_justchr.bed"
    line_bed=$(wc -l $tmp_bed)
    echo -e "bed: ${line_bed}"
    line_bed_cut=$(echo ${line_bed} | cut -d ' ' -f1 )

    if (( $line_maf_cut > 2 )) && (( $line_bed_cut > 1)); then
        echo $chr_clean

        # # cnee lrt ####
        in_sh="slurm_run_phylop_tmp_lrt_con_de_v2.sh"
        sbatch ${in_sh} ${tmp_bed} ${chr_clean} ${tmp_maf}

        # # cnee lrt branch
        tmp_br="HLlicCas1-HLacaPus1,HLtriMol2,HLcalAnn5-HLfloFus1"
        in_sh="slurm_run_phylop_tmp_branch_lrt_con_cnee_de_v2.sh"
        sbatch ${in_sh} ${tmp_br} ${tmp_bed} ${chr_clean} ${tmp_maf}

    fi
}

export -f runPhyloPbyChr
cut -f2 ${trans_matrix} | xargs -n 1 -P 10 -I {} bash -c 'runPhyloPbyChr "$@"' _ {}

# bash 2_Loop_chr_SUBMIT.sh > submit.log 2>&1

# check:
# ls |grep "^phylop" | grep ".err" | xargs cat > concat_err.log
# ls |grep "^phylop" | grep ".err"  |grep -nrl "ERROR: no features fall in alignment"
# chrUn_NW_020110153v1 is the one

# ls |grep "^phylop" | grep ".err" | grep -nrl "ERROR: no node named"
# then look up .out
# chrUn_NW_020110163v1
# chrUn_NW_020110159v1
# chrUn_NW_020110153v1
# chrUn_NW_020110162v1
# chr31
# chr32
