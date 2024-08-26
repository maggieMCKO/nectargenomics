#!/bin/bash


# wd: /home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/03_phyloacc/2024ver/nectar_nodollo/great

module purge
module load rev/11.06 bedtools/2.29.1 GREAT/1.5  

# GREAT sets: 
export great_dir="/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/03_phyloacc/2024ver/nectar_nodollo/great/"
export greatsets=(basalPlusExtension_proteincoding ) 

# Input ACC
# all cnees
export cneeAna_DIR="$HOME/Nectar/analysis/CNEEanalysis_2021/"
export acc_bed="${cneeAna_DIR}01b_cnee_prep/bed_outputs/galGal6_final_conserved_CNEEs.bed"

overlap_dir=overlap
mkdir -p ${overlap_dir}

## RUN START

# Loop through greatsets
for greatset in "${greatsets[@]}"
do
    echo "Processing greatset: $greatset"

    # full path of great regdom
    great_bed="${great_dir}galGal6_gene_tss_great_${greatset}.bed"

    output_file="${overlap_dir}/cnee_gene_assignment_${greatset}_allcnees.bed"

    bedtools intersect -wa -wb -a $acc_bed -b $great_bed > ${output_file}

done