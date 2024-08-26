#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --ntasks=24
#SBATCH --mem-per-cpu=5g
#SBATCH --partition=medium
#SBATCH -o log3_go_binomial_NectarSepcific_%A.o
#SBATCH -e log3_go_binomial_NectarSepcific_%A.e
#SBATCH -J go_binomial_NectarSepcific


# wd: /home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/03_phyloacc/2024ver/nectar_nodollo/great

module purge
# module load rev/11.06 bedtools/2.29.1 GREAT/1.5 
module load rev/20.12 rev/21.12 bedtools/2.29.1 great/1.5

# GO bed path
export go_bed_dir="/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/03_phyloacc/2024ver/nectar_nodollo/great/go/"
export onts=(biological_process molecular_function cellular_component) 

# GREAT sets: 
export greatsets=(basalPlusExtension_proteincoding ) # twoClosest_proteincoding

# Set
export acc_critera=(Loose NarrowBF2_gt3)

# Input ACC
export output_dir=/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/03_phyloacc/2024ver/nectar_nodollo/phyloacc-out-02-07-2024.11-55-54/results/acc_cnee_sets/
export acc_sets=(union congvergent2way congvergent3way congvergent4way honeyeaters hummingbirds parrots sunbirds congvergent_gteq2way congvergent_genelv_2way)

# Genome non-gap
export antigap_bed=galGal6_antigap.bed 

operation() {
    term=$1
    greatset=$2
    acc_set=$3
    num_cnee=$4
    ont=$5
    dir_tmp=$6

    # echo $greatset $acc_set $term $ont

    # term cleaned
	term_basename=$(basename "$term")
	# remove file extension
	term_cleaned="${term_basename%.*}"

	# full path of acc_set
	acc_bed="${output_dir}${acc_cri}_${acc_set}.bed" 

    numRegionsHit=$(bedtools intersect -u -wa -a $acc_bed -b $term | wc -l )

    # calculateBinomialP $regdoms_in $antigap_bed $numTotalRegions $numRegionsHit
	p=$(calculateBinomialP $term $antigap_bed $num_cnee $numRegionsHit)
	
	# great set | acc_critera | acc set | num of cnees in acc set | ont | term | number of acc. cnees overlap with the genes in the term | BinomialP
    echo -e "$greatset\t$acc_cri\t$acc_set\t$num_cnee\t$ont\t$term_cleaned\t$numRegionsHit\t$p" > "${dir_tmp}/${term_cleaned}.tmp"

}

export -f operation

## RUN START

# Loop through acc_critera
for acc_cri in "${acc_critera[@]}"
do

    # Loop through greatsets
    for greatset in "${greatsets[@]}"
    do
        echo "Processing greatset: $greatset"

        # Loop through acc_sets
        for acc_set in "${acc_sets[@]}"
        do
            # full path of acc_set
            acc_bed="${output_dir}${acc_cri}_${acc_set}.bed" 

            # Count the number of lines in the acc_set file
            num_cnee=$(wc -l < "$acc_bed")
            echo "  Acc set: $acc_set with num_cnee: $num_cnee"

            # Loop through onts
            for ont in "${onts[@]}"
            do
                # echo "    Ont: $ont"

                # Define bed directory for current greatset and ont
                bed_dir="${go_bed_dir}${greatset}/${ont}/bed/"

                # Define temporary directory for this run
                dir_tmp="tmp_go_${acc_cri}_${ont}_${greatset}_${acc_set}"
                mkdir -p "$dir_tmp"

                # Loop through terms (files in bed_dir)
                export greatset
                export acc_set
                export num_cnee
                export ont
                export dir_tmp
                export acc_cri
                find "$bed_dir" -type f -name "*.bed" | parallel -j+0 operation {} "$greatset" "$acc_set" "$num_cnee" "$ont" "$dir_tmp" "$acc_cri"

                # Concatenate all tmp files to the final output file
                mkdir -p great_BinomialP_go
                output_file="great_BinomialP_go/great_BinomialP_go_${ont}_${greatset}_${acc_cri}_${acc_set}.tsv"
                echo "output: ${output_file}"
                cat "${dir_tmp}"/*.tmp > "$output_file"

                # Remove the temporary directory and files
                rm -r "${dir_tmp}"
            done
        done
    done
done