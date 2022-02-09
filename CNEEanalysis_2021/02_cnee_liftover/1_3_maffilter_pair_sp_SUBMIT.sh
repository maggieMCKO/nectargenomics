#!/bin/bash

source "0_setup.sh"
mkdir -p d1_optionfiles d1_slurmscripts d1_logS
mkdir -p liftover_proc/0_maf/

export maffilter_optionfile="1_1_optionfile_tmp.maffilter"
export maffilter_sp_slurm="1_2_maffilter_pair_sp.sh"
export tmp_slurm_file=${maffilter_sp_slurm}
export tmp_maffilter_optionfile=${maffilter_optionfile}

MakeMafPair (){

	export sp=$1
    echo -e "current tmp_sp: ${sp}"

    # make slurm job script per sp
    tmp_slurm_file_stem=${tmp_slurm_file%%.sh}
    tmp_slurm_file_stem_sp="d1_slurmscripts/${tmp_slurm_file_stem}_${sp}.sh"
    echo -e "tmp_slurm_file_stem_sp: ${tmp_slurm_file_stem_sp}"
    sed "s/tmp/$sp/g" ${tmp_slurm_file} > ${tmp_slurm_file_stem_sp}

    # make optionfile per sp
    tmp_maffilter_optionfile_sp="d1_optionfiles/${tmp_maffilter_optionfile/tmp/$sp}"
    echo -e "tmp_maffilter_optionfile_sp: ${tmp_maffilter_optionfile_sp}"
    sed "s/VAR/$sp/g" ${tmp_maffilter_optionfile} > ${tmp_maffilter_optionfile_sp}

    # submit the job
    sbatch ${tmp_slurm_file_stem_sp}
    sleep 1
}

export -f MakeMafPair
cut -f1 ${specieslist} | grep "galGal6" | parallel MakeMafPair {}
# | grep -v "galGal6"
# bash 1_3_maffilter_pair_sp_SUBMIT.sh &> runninglog


# check progress
# find d1_logS  | grep ".o$" |grep -nrl "Elapsed" | xargs cat > mafFilter_prog.log



# check:
# find phylop_logs |grep "^phylop" | grep ".err" | xargs cat > concat_err.log
# grep -i error concat_err.log
# chr1 wasn't finished in 45 mins, used 66 mins

# ls |grep "^phylop" | grep ".err"  |grep -nrl "ERROR: no features fall in alignment"
