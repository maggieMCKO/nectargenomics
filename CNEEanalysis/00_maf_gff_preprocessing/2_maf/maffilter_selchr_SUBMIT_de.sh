#!/bin/bash

export cneeAna_DIR="/home/mpg08/mko/Nectar/analysis/CNEEanalysis/"
export maf_DIR="${cneeAna_DIR}CNEEanalysis/00_maf_gff_preprocessing/2_maf/"
export maffilter_optionfile="${maf_DIR}optionfile_tmp.maffilter"
export maffilter_selchr_slurm="${maf_DIR}maffilter_selchr_slurm_de.sh"

WD="${cneeAna_DIR}01_phylofit/"
trans_matrix=${WD}galgal6a_2315.6_acc3.tsv
Processed_DIR="${maf_DIR}Processed"
mkdir -p ${Processed_DIR}
log_DIR="${maf_DIR}logS"
mkdir -p ${log_DIR}

tmp_slurm_file=${maffilter_selchr_slurm}
tmp_maffilter_optionfile=${maffilter_optionfile}

chrs=($(cut -f2 ${trans_matrix}))

for chr in "${chrs[@]}"
do
    echo -e "current chr: ${chr}"

    # make slurm job script per chr
    tmp_slurm_file_stem=${tmp_slurm_file%%.sh}
    tmp_slurm_file_stem_chr="${tmp_slurm_file_stem}_${chr}.sh"
    echo -e "tmp_slurm_file_stem_chr: ${tmp_slurm_file_stem_chr}"
    sed "s/tmp/$chr/g" ${tmp_slurm_file} | sed "s/MF/MF_$chr/g" > ${tmp_slurm_file_stem_chr}

    # make optionfile per chr
    tmp_maffilter_optionfile_chr=${tmp_maffilter_optionfile/tmp/$chr}
    echo -e "tmp_maffilter_optionfile_chr: ${tmp_maffilter_optionfile_chr}"
    sed "s/tmp/$chr/g" ${tmp_maffilter_optionfile} > ${tmp_maffilter_optionfile_chr}

    # submit the job
    sbatch ${tmp_slurm_file_stem_chr}
    sleep 1

done
