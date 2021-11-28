#!/bin/bash
#SBATCH --time=0-01:30:00 # 1h for nectar final
#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=4g
#SBATCH --partition=medium
#SBATCH -o log1_4dsites_%A.out
#SBATCH -e log1_4dsites_%A.err
#SBATCH -J 4dsites
#SBATCH --mail-type=ALL        # BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mcko@orn.mpg.de
#SBATCH -C "scratch2"

##### SETUP INPUT and OUTPUT
export cneeAna_DIR="$HOME/Nectar/analysis/CNEEanalysis_2021/"
export input_DIR="${cneeAna_DIR}00_inputs/"
export input_tree="${input_DIR}NectarTree02.nw"
export input_maf="${input_DIR}multiz45way.maf.gz"

export maf_processed_DIR="${cneeAna_DIR}01a_phylofit/02_maf_prep/Processed/"

export gff_DIR="${cneeAna_DIR}01a_phylofit/01_gff_prep/"
export trans_matrix="${gff_DIR}galgal6a_2315.6_acc_chr_matrix.tsv" # (chr): for maf
export gff_processed_DIR="${gff_DIR}/galGal6_gff_split/" # (galGal6.chr): for gff # MUST

export WD="${cneeAna_DIR}01a_phylofit/03_run_phylofit/"
export tmp_DIR="${WD}01_msa/"
mkdir -p ${tmp_DIR}


export MYSCRATCH="/scratch/users/$USER/phylofit/"
mkdir -p ${MYSCRATCH}


## get species list
# module load ucsc/20160601
# gunzip -c ${input_maf} > "${input_DIR}multiz45way.maf"
# mafSpeciesList "${input_DIR}multiz45way.maf" "${input_DIR}sp_list.txt"

export sp_list=$(sed -z "s/\n/,/g;s/,$/\n/" "${input_DIR}sp_list.txt")
echo ${sp_list}

export concat_4ds="${WD}conca_4dSites.ss"
export neutral_model="${WD}nonconserved_4d"
export neutral_model_named="${WD}nonconserved_4d_named.mod"

##### SETUP programs
export phastPATH="$HOME/Tools/phast/bin/"
export phyloFit_PATH="${phastPATH}phyloFit"
export msa_view_PATH="${phastPATH}msa_view"
export tree_doctor_PATH="${phastPATH}tree_doctor"

##### RUN
### 1. extract 4d codons and 4d sites from an alignment (by chromosome)
Get4dCodonSitebyChr (){
    Chr=$1
    chr_clean=$(echo $Chr | sed 's/galGal6.//')

    # input maf
    chr_maf="${chr_clean}.maf"
    tmp_maf=$(find ${maf_processed_DIR} | grep ${chr_maf})
    line_maf=$(wc -l $tmp_maf)
    echo -e "maf: ${line_maf}"
    line_maf_cut=$(echo ${line_maf} | cut -d ' ' -f1 )

    # input gff
    # chr_gff="galGal6_fixed_${Chr}.gff"
    chr_gff="${Chr}.gff"
    tmp_gff=$(find ${gff_processed_DIR} | grep ${chr_gff})
    line_gff=$(wc -l $tmp_gff)
    echo -e "gff: ${line_gff}"
    line_gff_cut=$(echo ${line_gff} | cut -d ' ' -f1 )

    if (( $line_maf_cut > 2 )) && (( $line_gff_cut > 1)); then
        echo $chr_clean

        # output 4d codons and sites
        codons_4d="${out_DIR}${chr_clean}_4dCodons.ss"
        sites_4d="${out_DIR}${chr_clean}_4dSites.ss"

        # 4d codon
        ${msa_view_PATH} ${tmp_maf} --4d --features ${tmp_gff} > ${codons_4d}

        # 4d sites
        ${msa_view_PATH} ${codons_4d} --in-format SS --out-format SS --tuple-size 1 > ${sites_4d}

    fi
}

export -f Get4dCodonSitebyChr
cut -f2 ${trans_matrix} | parallel --memfree 4G --tmpdir ${MYSCRATCH} Get4dCodonSitebyChr {}

### 2. concatenate 4d sites
export non_zeros="${out_DIR}non_zeros"
mkdir -p ${non_zeros}

MoveNonZeros () {
    tmp_file=$1
    if [ -s ${tmp_file} ]
    then
         cp ${tmp_file} ${non_zeros}
    else
         echo "${tmp_file} empty"
    fi
}
export -f MoveNonZeros
find ${out_DIR} -maxdepth 1 | grep "4dSites.ss" | parallel MoveNonZeros {}

cd ${non_zeros}
${msa_view_PATH} --aggregate ${sp_list} --in-format SS --out-format SS --unordered-ss chr*_4dSites.ss > ${concat_4ds}
