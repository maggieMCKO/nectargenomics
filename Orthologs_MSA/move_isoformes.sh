#!/usr/bin/env bash
#

iso_file=$1
msa_dir=$2
msa_suffix=$3
out_dir=$4


while IFS= read -r line
do
	#echo $gene
	gene=$(awk '{print $1}' <(echo $line))
	transcript=$(awk '{print $2}' <(echo $line))
	msa_file=$msa_dir/$transcript.$msa_suffix
	echo $msa_file
	if [ -s $msa_file ]; then
		
		cp $msa_dir/$transcript.$msa_suffix $out_dir/$gene/
	fi
done < "$iso_file"
