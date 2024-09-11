#!/bin/bash


# wd: /home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/03_phyloacc/2024ver/nectar_nodollo/great

module purge
# module load rev/11.06 bedtools/2.29.1 GREAT/1.5 
module load rev/20.12 rev/21.12 bedtools/2.29.1 great/1.5

# GREAT sets: 
export great_dir="/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/03_phyloacc/2024ver/nectar_nodollo/great/"
export greatsets=(basalPlusExtension_proteincoding ) # /home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/03_phyloacc/oldway_rerun_rmPar/post_phyloacc/great

# Set
export acc_critera=(Loose NarrowBF2_gt3)

# Input ACC
export output_dir=/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/03_phyloacc/2024ver/nectar_nodollo/phyloacc-out-02-07-2024.11-55-54/results/acc_cnee_sets_v7/
export acc_sets=(union congvergent2way congvergent3way congvergent4way honeyeaters hummingbirds parrots sunbirds congvergent_gteq2way congvergent_gteq3way)

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

            output_file="${overlap_dir}/cnee_gene_assignment_${greatset}_${acc_cri}_${acc_set}_v72.bed"

            bedtools intersect -wa -wb -a $acc_bed -b $great_bed > ${output_file}

        done
    done
done