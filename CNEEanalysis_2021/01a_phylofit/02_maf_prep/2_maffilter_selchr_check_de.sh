#!/bin/bash

export cneeAna_DIR="$HOME/Nectar/analysis/CNEEanalysis_2021/"
export WD="${cneeAna_DIR}01a_phylofit/02_maf_prep/"
export log_DIR="${WD}logS"



out_logs=($(find ${log_DIR}| grep "\.o" ))          # array

for log in "${out_logs[@]}"
do
    if grep -q "MafFilter's done. Bye." ${log}; then # -z: if string is null
        : # do nothing
    else
        echo -e "${log} did not finished"
    fi

done

# bash 2_maffilter_selchr_check_de.sh &> check_progress.txt
