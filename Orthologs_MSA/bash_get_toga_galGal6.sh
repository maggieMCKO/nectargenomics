#!/usr/bin/bash

ref=galGal6
#defQuery=HLacaPus1
trans=$1
db=HLstrHab1
projections=/projects/project-osipova/NectarivoryProject/TOGA_ref_species_ncbi/TOGA_galGal6_new_run/${db}_toga_new_run/projections.fasta
ortho_file=/projects/project-osipova/NectarivoryProject/TOGA_ref_species_ncbi/TOGA_galGal6_new_run/${db}_toga_new_run/orthology_type.tsv
out_dir=orthologs_fa_from_toga/
out_fasta=$out_dir/$trans/$trans.$ref.fa

## prepare output dir structure
mkdir -p $out_dir
mkdir -p $out_dir/$trans


## check if all required files are there
if [ ! -f $projections ]; then
	echo "projections file $projections was not found!"
        exit 1
fi

if [ ! -f $ortho_file ]; then
        echo "orthology type file $ortho_file was not found!"
	exit 1
fi

## get query transcript name
#queryTrans=$(grep $trans $ortho_file | cut -f4)
#queryTrans=${queryTrans#GREY_}
#echo $queryTrans
#echo $db

## extract sequences of this transcript from toga projections
get_seq_by_name_fasta.py -f $projections  -n ref_${trans} -p "${ref};" > $out_fasta
