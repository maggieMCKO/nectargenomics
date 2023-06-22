#!/usr/bin/env bash
#

gene_file=$1
msa_final_dir=$2
absrel_dir=$3
# only these transcripts are to be considered
qual_trans=$4
# suffix to remove
suff=$5

for gene in $(cut -f1 $gene_file);
do
	if [ -d $msa_final_dir/$gene/  ]; then

		for t in $(ls $msa_final_dir/$gene/)
		do
			if  grep -q "${t%$suff}" $qual_trans ; then
				ans_tree=${t%$suff}.ans.tree
				out_file="out.${t%$suff}.tsv"
				echo "run_phyphy.py -a $msa_final_dir/$gene/$t -t $absrel_dir/$gene/${t%$suff}/${t%$suff}.ans.tree -f -m absrel -o $absrel_dir/$gene/${t%$suff}/$out_file"
			fi
		done
	fi
done

