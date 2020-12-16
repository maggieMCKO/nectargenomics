#!/usr/bin/env bash
#
# prune clades other than hummingbirds!

gene_file=$1
msa_final_dir=$2
absrel_dir=$3
tree=$4
# only these transcripts have to be considered
qual_trans=$5
# suffix to remove
suff=$6


for gene in $(cut -f1 $gene_file);
do
	if [ -d $msa_final_dir/$gene/  ]; then
		mkdir -p $absrel_dir/$gene

		for t in $(ls $msa_final_dir/$gene/)
		do
			mkdir -p $absrel_dir/$gene/${t%$suff}
			
			if  grep -q "${t%$suff}" $qual_trans ; then
				# prune tree and add ancestral branches
				echo "label.py <( tree_doctor -a -n -P \$(grep '>' $msa_final_dir/$gene/$t | sed 's/>//g' |  grep -v 'HLlicCas1\|HLtriMol2\|HLphyNov1' | tr '\n' ',') $tree | sed 's/-/_/g') HLcalAnn5_HLfloFus1 > $absrel_dir/$gene/${t%$suff}/${t%$suff}.ans.tree"
			fi
		
		done
	fi
done

