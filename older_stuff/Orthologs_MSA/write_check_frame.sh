#!/usr/bin/env bash
#

trans=$1
fasta=$2
frame_symbol=$3


for fa_name in $(grep '>' $fasta);
do
    name=${fa_name#>}
    if  [ $(check_frame_fasta.py -f $fasta -s $name -n $frame_symbol) ]; then
        frame=$frame_symbol
    else
        frame="0"
    fi

    echo -e "$trans\t$name\t$frame"
done
