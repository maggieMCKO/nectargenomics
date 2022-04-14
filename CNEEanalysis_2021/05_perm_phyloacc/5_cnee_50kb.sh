#!/bin/bash




# WD: "/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/05_perm_phyloacc/1_prepare_new/"

module purge
module load bedtools/2.29.1 anaconda2/2019.10
source activate R_phylo


mkdir -p postPhyloAcc
cd postPhyloAcc/

mkdir linear

# set 50 kb window around genes
bedtools slop -i galGal.ncbigenes.sorted.bed -g galGal6.chrom.sizes -b 50000 > linear/galGal.ncbislop_50kb.bed



bedtools intersect -a "accelerated_cnees_noheader.bed" -b "linear/galGal.ncbislop_50kb.bed" -wa -wb > linear/cnee_ncbigene50kb_acc_intersect.bed
# awk '{print NF}' linear/cnee_ncbigene50kb_acc_intersect.bed  | sort -u
bedtools groupby -i linear/cnee_ncbigene50kb_acc_intersect.bed -g 10,11,12,13 -c 7 -o max > linear/cnee_ncbigene50kb_acc_maxExp.bed