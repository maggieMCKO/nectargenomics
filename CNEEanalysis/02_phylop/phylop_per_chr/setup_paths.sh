#!/bin/bash
##### SETUP INPUT and OUTPUT
export cneeAna_DIR="/home/mpg08/mko/Nectar/analysis/CNEEanalysis/"
export input_tree="${cneeAna_DIR}NectarTree01.nw"
export chicken_anno="${cneeAna_DIR}GCF_000002315.6_GRCg6a_genomic_chrtabFixed2.gff" # chr need to match maf's chr with species (galgal.chrx)

export maf_DIR="${cneeAna_DIR}MAF/"
export input_maf="${maf_DIR}multi.anno.maf"
export maf_processed_DIR="${maf_DIR}Processed"

export phylofit_DIR="${cneeAna_DIR}01_phylofit/01_msa3/"
export neutral_model="${phylofit_DIR}nonconserved_4d.mod"

export phylop_DIR="${cneeAna_DIR}02_phylop/"
export features_bed="${phylop_DIR}galGal6_final_merged_CNEEs_named_fixchr_justchr.bed"
# features should be just chr, but maf species.chr
# features_out="${phylop_DIR}features.out"

export WD="${cneeAna_DIR}02_phylop/check_detail/"
export neutral_model_named="${WD}nonconserved_4d_named.mod"

export trans_matrix="${cneeAna_DIR}01_phylofit/galgal6a_2315.6_acc3.tsv"  #galgal6.chrx
export trans_matrix2="${cneeAna_DIR}01_phylofit/galgal6a_2315.6_acc5.tsv"
# just chr for bed

export CNEEbyChr="${phylop_DIR}CNEE_byChr/"

##### SETUP programs
export phastPATH="/home/mpg08/mko/Tools/phast/bin/"
export phyloFit_PATH="${phastPATH}phyloFit"
export msa_view_PATH="${phastPATH}msa_view"
export phyloP_PATH="${phastPATH}phyloP"
export tree_doctor_PATH="${phastPATH}tree_doctor"
