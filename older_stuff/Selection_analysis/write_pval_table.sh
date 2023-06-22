#!/usr/bin/env bash
#

iso_file=$1
absrel_dir=$2


while IFS= read -r line
do
	gene=$(echo $line | cut -d' ' -f1)
	trans=$(echo $line | cut -d' ' -f2)
	if [ -s "$absrel_dir/$gene/$trans/pval.$trans.tsv"  ]; then
		sed "s/$/\t$trans\t$gene/" "$absrel_dir/$gene/$trans/pval.$trans.tsv"
	fi
done < $iso_file
