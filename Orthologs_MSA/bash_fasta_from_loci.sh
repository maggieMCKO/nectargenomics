#!/usr/bin/env bash
#

## initialize params
trans=$1
db=$2
projections_bed=$3
anno_bed=cds_annotations/cds.$db.anno.bed
out_dir=orthologs_fasta_by_transcript/
out_fasta=$out_dir/$trans/$trans.$db.dna.fa
fraction=0.5


## prepare output dir structure
mkdir -p $out_dir
mkdir -p $out_dir/$trans


## get TOGA orthologous loci for transcripts
#echo "Getting TOGA orthologous loci for transcript $trans ...."

if [ -f $projections_bed  ]; then
	toga_loci=$(grep $trans $projections_bed) 
else
	echo "projections file $projections_bed was not found!"
	exit 1
fi


if [ "$toga_loci" == "" ]; then
	echo "didn't find anything for $trans in $projections_bed file!"
	exit 0
fi


## get corresponding transcript structure from annotation file 
#echo "Getting corresponding transcript structure from annotation file ...."

if [ -f $anno_bed ]; then
	#echo "running: bedtools intersect -split -a $anno_bed -b <(echo "\$toga_loci") -wo -f $fraction"
	anno_loci=$(bedtools intersect -split -a $anno_bed -b <(echo "$toga_loci") -wo -f $fraction)
	#echo "$anno_loci"
	#exit 0
	if [ "$anno_loci" == ""  ]; then
		echo "there's nothing for $toga_loci in $anno_bed!"
		exit 0
	else
		best_anno_name=$(python /projects/project-osipova/PythonScripts/get_best_matching_transcript.py -t <(echo "$anno_loci"))
		best_anno_trans=$(grep -w $best_anno_name $anno_bed)
	fi
else
	echo "annotaiton file $anno_bed was not found!"
        exit 1
fi


## get fasta sequence for this transcript
#echo "Getting DNA fasta for the transcript ...."
#echo "$best_anno_name"
#echo "$best_anno_trans"
bedtools getfasta -name -split -s -fi /projects/hillerlab/genome/gbdb-HL/$db/$db.fa -bed <(echo "$best_anno_trans") -fo $out_fasta

#bedToGenePred <(echo $best_anno_trans) stdout | genePredToProt -cdsFa=$out_fasta stdin /projects/hillerlab/genome/gbdb-HL/$db/$db.2bit /dev/null

