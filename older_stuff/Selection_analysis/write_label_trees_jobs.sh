#!/usr/bin/env bash
#

gene_file=$1
msa_final_dir=$2
absrel_dir=$3
fg_branches=$4
old_absrel_dir=$5
# only these transcripts have to be considered
qual_trans=$6
# suffix to remove
suff=$7


for gene in $(cut -f1 $gene_file);
do
	if [ -d $msa_final_dir/$gene/  ]; then
		mkdir -p $absrel_dir/$gene

		for t in $(ls $msa_final_dir/$gene/)
		do
			mkdir -p $absrel_dir/$gene/${t%$suff}
			
			if  grep -q "${t%$suff}" $qual_trans ; then
				# label Foreground branches in the existing tree
				echo "label.py  $old_absrel_dir/$gene/${t%$suff}/${t%$suff}.ans.tree $fg_branches > $absrel_dir/$gene/${t%$suff}/${t%$suff}.ans.tree"
			fi
		
		done
	fi
done

