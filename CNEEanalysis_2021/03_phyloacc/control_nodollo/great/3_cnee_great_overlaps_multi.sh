#!/bin/bash
# wd: /home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/03_phyloacc/2024ver/control_nodollo/great

module purge
module load rev/20.12 rev/21.12 bedtools/2.29.1 great/1.5

# GREAT sets: 
export great_dir="/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/03_phyloacc/2024ver/nectar_nodollo/great/"
export greatsets=(basalPlusExtension_proteincoding ) 

# Set
export acc_critera=( ControlLoose )

# Input ACC
export output_dir=/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/03_phyloacc/2024ver/control_nodollo/phyloacc-out-02-20-2024.01-26-12/results/acc_cnee_sets/
export acc_sets=(union congvergent2way congvergent3way congvergent4way swifts falcons lyrebirds passerides congvergent_gteq2way congvergent_genelv_gteq2way)

## RUN START

# Loop through acc_critera
for acc_cri in "${acc_critera[@]}"
do
    overlap_dir=overlap/${acc_cri}
    mkdir -p ${overlap_dir}

    # Loop through greatsets
    for greatset in "${greatsets[@]}"
    do
        echo "Processing greatset: $greatset"

        # Loop through acc_sets
        for acc_set in "${acc_sets[@]}"
        do
            # full path of acc_set
            acc_bed="${output_dir}${acc_cri}_${acc_set}.bed" 

            # full path of great regdom
            great_bed="${great_dir}galGal6_gene_tss_great_${greatset}.bed"

            output_file="${overlap_dir}/cnee_gene_assignment_${greatset}_${acc_cri}_${acc_set}.bed"

            bedtools intersect -wa -wb -a $acc_bed -b $great_bed > ${output_file}

        done
    done
done