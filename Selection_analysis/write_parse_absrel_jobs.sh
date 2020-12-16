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
				out_json="out.${t%$suff}.tsv"
				pval_tab="pval.${t%$suff}.tsv"	
				echo "parse_absrel_json.py -a 'Corrected P-value' -j $absrel_dir/$gene/${t%$suff}/$out_json > $absrel_dir/$gene/${t%$suff}/$pval_tab"
			fi
		done
	fi
done

