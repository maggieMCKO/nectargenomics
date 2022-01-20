#!/bin/bash
#SBATCH --time=00:30:00 # run time: 10m
#SBATCH --ntasks=10
#SBATCH --mem=5g
#SBATCH --partition=medium
#SBATCH -o submitphylop_%A.out
#SBATCH -e submitphylop_%A.err
#SBATCH -J submitphylop

export cneeAna_DIR="../"
export trans_matrix_for_maf="${cneeAna_DIR}01a_phylofit/01_gff_prep/galgal6a_2315.6_acc_chr_matrix.tsv" # (chr)

export maf_processed_DIR="${cneeAna_DIR}01a_phylofit/02_maf_prep/Processed/"


# 1. run phyloP by chromosome
runPhyloPbyChr (){
    Chr=$1

    # input maf
    chr_maf="${Chr}.maf"
    tmp_maf=$(find ${maf_processed_DIR} | grep ${chr_maf})
    line_maf=$(wc -l $tmp_maf)
    echo -e "maf: ${line_maf}"
    line_maf_cut=$(echo ${line_maf} | cut -d ' ' -f1 )

    # input bed
    tmp_bed="cnee_byChr/galGal6_${Chr}_justchr.bed"
    line_bed=$(wc -l $tmp_bed)
    echo -e "bed: ${line_bed}"
    line_bed_cut=$(echo ${line_bed} | cut -d ' ' -f1 )

    if (( $line_maf_cut > 2 )) && (( $line_bed_cut > 1)); then
        echo $Chr

        # cnee lrt con ####
        # in_sh="4_2_phylop_tmplet_lrt_con.sh"
        in_sh=$2
        echo -e "in_sh: $in_sh"
        sbatch ${in_sh} ${tmp_bed} ${Chr} ${tmp_maf}

        sleep 1
    fi
}

export -f runPhyloPbyChr

# chr1: becasue chr1 needs more run time
runPhyloPbyChr "chr1" "4_2_phylop_tmplet_lrt_con_chr1.sh"
# other chr
cut -f2 ${trans_matrix_for_maf} |grep -v "chr1$" | parallel runPhyloPbyChr {} 4_2_phylop_tmplet_lrt_con.sh



# check:
# find phylop_logs |grep "^phylop" | grep ".err" | xargs cat > concat_err.log
# grep -i error concat_err.log
# chr1 wasn't finished in 45 mins, used 66 mins

# ls |grep "^phylop" | grep ".err"  |grep -nrl "ERROR: no features fall in alignment"
# then look up .out
# phylop_logs/phylop_10927043.err chrUn_NW_020110147v1 v
# phylop_logs/phylop_10927042.err chrUn_NW_020110134v1 v
# phylop_logs/phylop_10927058.err chrUn_NW_020110156v1 v
# phylop_logs/phylop_10927044.err chrUn_NW_020110150v1 v
# phylop_logs/phylop_10927051.err chrUn_NW_020109947v1 v

# concatenate phyloP_output
# cat phyloP_outputs/*.tsv > phylop_allChr_lrt_con.tsv
