#!/bin/bash
#SBATCH --time=0-06:00:00
#SBATCH -n 1 # or --ntasks=1
#SBATCH --mem=20g
#SBATCH --partition=medium
#SBATCH -o phylopSC_%A.out
#SBATCH -e phylopSC_%A.err
#SBATCH -J phylopSC
#SBATCH --mail-type=ALL        # BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mcko@orn.mpg.de

##### SETUP INPUT and OUTPUT
cneeAna_DIR="/home/mpg08/mko/Nectar/analysis/CNEEanalysis/"
nput_tree="${cneeAna_DIR}NectarTree01.nw"
chicken_anno="${cneeAna_DIR}GCF_000002315.6_GRCg6a_genomic_chrtabFixed2.gff" # chr need to match maf's chr with species (galgal.chrx)

maf_DIR="${cneeAna_DIR}MAF/"
input_maf="${maf_DIR}multi.anno.maf"
maf_processed_DIR="${maf_DIR}Processed"

phylofit_DIR="${cneeAna_DIR}01_phylofit/01_msa3/"
neutral_model="${phylofit_DIR}nonconserved_4d.mod"

phylop_DIR="${cneeAna_DIR}02_phylop/"
features_bed="${phylop_DIR}galGal6_final_merged_CNEEs_named_fixchr2.bed"
features_out="${phylop_DIR}features.out"

WD="${cneeAna_DIR}02_phylop/subtree2/"
neutral_model_named="${WD}nonconserved_4d_named.mod"

tmp_subtree="HLcalAnn5-HLfloFus1"
echo ${tmp_subtree}

##### SETUP programs
export phastPATH="/home/mpg08/mko/Tools/phast/bin/"
export phyloFit_PATH="${phastPATH}phyloFit"
export msa_view_PATH="${phastPATH}msa_view"
export phyloP_PATH="${phastPATH}phyloP"
export tree_doctor_PATH="${phastPATH}tree_doctor"

##### RUN
### Run phyloP with the branch options.
${phyloP_PATH} --method LRT --subtree ${tmp_subtree} --mode CONACC --features ${features_bed} --msa-format MAF ${neutral_model_named} ${input_maf} > "subtree_${tmp_subtree}_basebybase_cnee"
