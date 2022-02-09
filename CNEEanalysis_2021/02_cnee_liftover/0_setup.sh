#!/bin/bash
module load ucsc/20160601
module load bedtools/2.29.1
# ls /usr/product/bioinfo/SL_7.0/BIOINFORMATICS/UCSC/20160601 | grep -i psl
export ucsc_path=$HOME/Tools/ucsc_liftover/bin/
# export seqtk_path=$HOME/Tools/seqtk/seqtk # not necessary

export cneeAna_DIR="../"
export Ref_species="galGal6" # yes: 1

# chr list
# export trans_matrix_for_maf="${cneeAna_DIR}01a_phylofit/01_gff_prep/galgal6a_2315.6_acc_chr_matrix.tsv" # (chr) # yes: 1

# chicken anno
# export fixed_chicken_anno="${cneeAna_DIR}00_inputs/GCF_000002315.6_GRCg6a_genomic_chrtabFixed2.gff" # yes: 1

# maf per chr
# export maf_processed_DIR="${cneeAna_DIR}01a_phylofit/02_maf_prep/Processed/" # yes: 1

# species list
export specieslist="${cneeAna_DIR}00_inputs/sp_list.txt" # yes: 1,2

# final CNEE bed
export filtered_final="../01b_cnee_prep/bed_outputs/galGal6_final_conserved_CNEEs.bed"
# yes: 1,2


export WD="./" # /home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/02_cnee_liftover
# export cneePSL="${WD}cnee.psl" # yes: 1
# export tmp_pslMappedTogalGal_dir="${WD}d5_CneeMappedToGg6_PSL_DIR/" # yes: 1,2


# export Renamed_DIR="${cneeAna_DIR}00_inputs/genomes/" # no: 1; yes: 2

# 3
export scratch_DIR="/scratch/users/$USER/phyloacc/"
mkdir -p ${scratch_DIR} # yes: 1

export aligned_DIR="${scratch_DIR}aligned/"
mkdir -p ${aligned_DIR} # no: 1

export unaligned_DIR="${scratch_DIR}unaligned/"
mkdir -p ${unaligned_DIR} # no: 1
