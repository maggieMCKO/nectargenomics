#!/bin/bash
#SBATCH --time=1-07:00:00 # used 557 mins
#SBATCH -N 4 # 4 for 96
#SBATCH --ntasks=96
#SBATCH --mem-per-cpu=4g
#SBATCH --partition=medium
#SBATCH -o log_csubst_t_%A.o
#SBATCH -e log_csubst_r_%A.e
#SBATCH -J csubst_r

# wd: /home/mpg08/mko/Nectar/analysis/csubst/01_run_csubst

export seqkit_PATH="$HOME/Tools/seqkit"

export wd="$HOME/Nectar/analysis/csubst/01_run_csubst/"
export input_dir="$HOME/Nectar/analysis/csubst/00_inputs/"
export nt_dir="${input_dir}one2one_extracted_toga/"
export gene_list="${input_dir}isoforms.one2ones.checked.tsv"
export tree_file="${input_dir}nonconserved_4d.nw"
export foreground="${input_dir}foreground3_nopar.tsv"
export scratch_DIR="/scratch/users/$USER/csubst/"


export proj="out_v3"

mkdir -p nt_align_renamed ${proj}
mkdir -p ${scratch_DIR}

runcsubst () {
	export tmp_transcipt=$1

	module load anaconda2/2019.10
	source activate env_csubst

	## nt align
	export nt_align=$(find ${nt_dir} | grep ${tmp_transcipt})
	# rename fa headers
	export tmp_intermediate_dir="${wd}nt_align_renamed/${tmp_transcipt}/"
	mkdir -p ${tmp_intermediate_dir}

	if [ ! -d ${tmp_intermediate_dir} ] ; then
			# if dir doesnt exist

			## rename fa headers
			export nt_renamed_align="${tmp_intermediate_dir}${tmp_transcipt}_renamed.fa"
			${seqkit_PATH} replace -p "\s.+" ${nt_align} | ${seqkit_PATH} replace -p "vs_" -r "" | ${seqkit_PATH} replace -p "REFERENCE" -r "galGal6" | ${seqkit_PATH} seq -w 0 > ${nt_renamed_align}

			## prune tree (remove tips not in the alignment)
			python prune.py ${tree_file} ${nt_renamed_align} ${tmp_intermediate_dir} ${tmp_transcipt}
			# export pruned_tree=$(find ${tmp_intermediate_dir} | grep "pruned.nw")
			export pruned_tree="${tmp_intermediate_dir}${tmp_transcipt}_pruned.nw"

	else
		if [ "$(ls -A ${tmp_intermediate_dir})" ] ; then
		# if dir is not empty

		export nt_renamed_align="${tmp_intermediate_dir}${tmp_transcipt}_renamed.fa"
		export pruned_tree="${tmp_intermediate_dir}${tmp_transcipt}_pruned.nw"

			if [ ! -f ${nt_renamed_align} ] ; then
				# if file doesnt exist
				## rename fa headers
				export nt_renamed_align="${tmp_intermediate_dir}${tmp_transcipt}_renamed.fa"
				${seqkit_PATH} replace -p "\s.+" ${nt_align} | ${seqkit_PATH} replace -p "vs_" -r "" | ${seqkit_PATH} replace -p "REFERENCE" -r "galGal6" | ${seqkit_PATH} seq -w 0 > ${nt_renamed_align}
			fi

			if [ ! -f ${pruned_tree} ] ; then
				# if file doesnt exist
				## prune tree (remove tips not in the alignment)
				python prune.py ${tree_file} ${nt_renamed_align} ${tmp_intermediate_dir} ${tmp_transcipt}
				# export pruned_tree=$(find ${tmp_intermediate_dir} | grep "pruned.nw")
				export pruned_tree="${tmp_intermediate_dir}${tmp_transcipt}_pruned.nw"
			fi
		fi
	fi

}
export -f runcsubst
cut -f1 ${gene_list} | parallel --max-procs ${SLURM_NTASKS} --memfree 4G runcsubst {}
