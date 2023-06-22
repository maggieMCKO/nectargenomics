#!/bin/bash

source setup_paths.sh
echo $features_bed

splitbedbyChr (){
    Chr=$1
    echo $Chr

    # input bed
    tmp_bed="${CNEEbyChr}galGal6_${Chr}_justchr.bed"
    awk '/'$Chr'\t/' ${features_bed} > ${tmp_bed}
}

export -f splitbedbyChr
cut -f2 ${trans_matrix2} | xargs -n 1 -P 10 -I {} bash -c 'splitbedbyChr "$@"' _ {}

# bash 1_splitCneeBedbyChr.sh > splitcnee.log 2>&1
# to check: ls | xargs wc -l
