#!/bin/bash
#SBATCH --time=0-4:00:00
#SBATCH -N 1
#SBATCH --ntasks=2
#SBATCH --mem=24g
#SBATCH --partition=fat
#SBATCH -o logs/log_csubst_31b_%A_%a.o # <------- change here
#SBATCH -e logs/log_csubst_31b_%A_%a.e # <------- change here
#SBATCH -J csubst
#SBATCH -C "scratch2"
#SBATCH --array 0-1968 # <------- change here


# wd: /home/mpg08/mko/Nectar/analysis/csubst/01_run_csubst

export seqkit_PATH="$HOME/Tools/seqkit"

export wd="$HOME/Nectar/analysis/csubst/01_run_csubst/"
export input_dir="$HOME/Nectar/analysis/csubst/00_inputs/"
export nt_dir="${input_dir}one2one_extracted_toga/"
export gene_list="${input_dir}isoforms.one2ones.checked.tsv"
export tree_file="${input_dir}nonconserved_4d.nw"
export foreground="${input_dir}foreground3.tsv" # run1 and run2
export foreground="${input_dir}foreground3_nopar.tsv"
export scratch_DIR="/scratch/users/$USER/csubst/"
export batch_DIR="${wd}batches13/" # <------- change here

export proj="out_v3"

mkdir -p nt_align_renamed ${proj} logs
mkdir -p ${scratch_DIR}

module load anaconda2/2019.10
source activate env_csubst

runcsubst () {
	export tmp_transcipt=$1
	echo -e "tmp_transcipt: ${tmp_transcipt}"

	## nt align
	export nt_align=$(find ${nt_dir} | grep ${tmp_transcipt})
	# rename fa headers
	export tmp_intermediate_dir="${wd}nt_align_renamed/${tmp_transcipt}/"
	mkdir -p ${tmp_intermediate_dir}

	echo "tmp_intermediate_dir: ${tmp_intermediate_dir}"

	if [ -d ${tmp_intermediate_dir} ] ; then
			# if dir exists
		echo -e "tmp_intermediate_dir exists"
		if [ "$(ls -A ${tmp_intermediate_dir})" ] ; then
			# if dir is not empty
			# echo "tmp_intermediate_dir not empty"
			export nt_renamed_align="${tmp_intermediate_dir}${tmp_transcipt}_renamed.fa"
			export pruned_tree="${tmp_intermediate_dir}${tmp_transcipt}_pruned.nw"
			echo "pruned_tree: ${pruned_tree}"
			if [ -f ${pruned_tree} ] ; then
				# if file exists
				echo -e "pruned_tree exists"

				# Run3
				export tmp_dir="${scratch_DIR}${proj}/${tmp_transcipt}/run3/"
				mkdir -p ${tmp_dir}
				export log="${tmp_dir}/log"

				export tmp_dir_wd="${wd}${proj}/${tmp_transcipt}/"
				mkdir -p ${tmp_dir_wd}

				if [ ! -d ${tmp_dir} ] ; then
					# if dir doesnt exist
					echo -e "run3: ${tmp_transcipt}"

					cd ${tmp_dir}
					csubst analyze \
					--alignment_file ${nt_renamed_align} \
			        --rooted_tree_file ${pruned_tree} \
			        --foreground ${foreground} \
			        --fg_exclude_wg no \
			        --max_arity 10 \
			        --exhaustive_until 2 &> ${log}

					cp -r "${tmp_dir}" "${tmp_dir_wd}"

				else
					if [ "$(ls -A ${tmp_dir})" ] ; then
						# if dir is not empty
						echo -e "tmp_dir (run3) ${tmp_dir} not empty"

						if ! grep -q "csubst analyze: Time elapsed" ${log}; then
							# doesnt found PATTERN

							echo -e "RERUN-1: run3: ${tmp_transcipt}"

							cd ${tmp_dir}
							csubst analyze \
							--alignment_file ${nt_renamed_align} \
					        --rooted_tree_file ${pruned_tree} \
					        --foreground ${foreground} \
					        --fg_exclude_wg no \
					        --max_arity 10 \
					        --exhaustive_until 2 &> ${log}

							cp -r "${tmp_dir}" "${tmp_dir_wd}"
						else
							cp -r "${tmp_dir}" "${tmp_dir_wd}" # add 0618
						fi
					else
					# if dir is empty
						echo -e "RERUN-2: run3: ${tmp_transcipt}"

						cd ${tmp_dir}
						csubst analyze \
						--alignment_file ${nt_renamed_align} \
				        --rooted_tree_file ${pruned_tree} \
				        --foreground ${foreground} \
				        --fg_exclude_wg no \
				        --max_arity 10 \
				        --exhaustive_until 2 &> ${log}

						cp -r "${tmp_dir}" "${tmp_dir_wd}"
					fi
				fi
			fi # pruned tree exists
		fi # nt_renamed not empty
	fi # nt_renamed dir exists

}
export -f runcsubst

printf -v BATCH "%04d" "${SLURM_ARRAY_TASK_ID}"
cut -f1 ${batch_DIR}${BATCH} | parallel  runcsubst {}
# echo -e "${BATCH} finished"
