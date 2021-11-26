#!/bin/bash

export cneeAna_DIR="$HOME/Nectar/analysis/CNEEanalysis_2021/"

export trans_matrix="${cneeAna_DIR}01a_phylofit/01_gff_prep/galgal6a_2315.6_acc_chr_matrix.tsv" # Must use acc5 format (chr)
export maf_DIR="${cneeAna_DIR}00_inputs/"
export maffilter_optionfile="1_1_optionfile_tmp.maffilter"
export maffilter_selchr_slurm="1_2_maffilter_selchr_slurm_de.sh"
export tmp_slurm_file=${maffilter_selchr_slurm}
export tmp_maffilter_optionfile=${maffilter_optionfile}

export WD="${cneeAna_DIR}01a_phylofit/02_maf_prep/"
mkdir -p ${WD}

export Processed_DIR="${WD}Processed"
mkdir -p ${Processed_DIR}
export log_DIR="${WD}logS"
mkdir -p ${log_DIR}
export tmp_script_DIR="${WD}running_scripts"
mkdir -p ${tmp_script_DIR}

prep (){

	export chr=$1
    echo -e "current chr: ${chr}"

    # make slurm job script per chr
    tmp_slurm_file_stem=${tmp_slurm_file%%.sh}
    tmp_slurm_file_stem_chr="${tmp_script_DIR}/${tmp_slurm_file_stem}_${chr}.sh"
    echo -e "tmp_slurm_file_stem_chr: ${tmp_slurm_file_stem_chr}"
    sed "s/tmp/$chr/g" ${tmp_slurm_file} | sed "s/MF/$chr/g" > ${tmp_slurm_file_stem_chr}

    # make optionfile per chr
    tmp_maffilter_optionfile_chr="${tmp_script_DIR}/${tmp_maffilter_optionfile/tmp/$chr}"
    echo -e "tmp_maffilter_optionfile_chr: ${tmp_maffilter_optionfile_chr}"
    sed "s/tmp/$chr/g" ${tmp_maffilter_optionfile} > ${tmp_maffilter_optionfile_chr}

    # submit the job
    sbatch ${tmp_slurm_file_stem_chr}
    sleep 1
}

export -f prep
cut -f2 ${trans_matrix} | parallel prep {}

# bash 1_3_maffilter_selchr_SUBMIT_de.sh &> runninglog
