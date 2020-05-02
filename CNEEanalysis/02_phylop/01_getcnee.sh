#!/bin/bash

#### Building a set of consensus CNEEs from literature ####
# based on: https://github.com/sjswuitchik/duck_comp_gen/blob/master/03a_cnee_analysis/01_assemble_cnees.sh

galgal6_gff=/home/mpg08/mko/Nectar/analysis/CNEEanalysis/GCF_000002315.6_GRCg6a_genomic.gff
replace_chrs_perl=/home/mpg08/mko/Nectar/analysis/CNEEanalysis/01_phylofit/replace_chrs.pl
trans_matrix=/home/mpg08/mko/Nectar/analysis/CNEEanalysis/01_phylofit/galgal6a_2315.6_acc5.tsv
final_cnee=galGal6_final_merged_CNEEs_named_fixchr.bed


# extract .gz
# ls | grep .gz| xargs -n 1 -P 4 gunzip

module purge
module load BEDTOOLS/2.29.1

# get ce lengths
awk '{print $3-$2, "\tLowe"}' galGal6_Lowe_merged.bed  >> ce.lengths
awk '{print $3-$2, "\tSackton"}' galGal6_Sackton_merged.bed  >> ce.lengths
awk '{print $3-$2, "\tUCSC"}' galGal6_UCSC_merged.bed  >> ce.lengths
awk '{print $3-$2, "\tCraig"}' galGal6_Craig_merged.bed  >> ce.lengths

cat galGal6_*_merged.bed | bedtools sort -i - | bedtools merge -i - | awk '{if ($3-$2 >= 50) print $0}' > galGal6_all_merged.bed

awk '{print $3-$2}' galGal6_all_merged.bed > galGal6_allce.lengths

# get exons
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.gff.gz
# gunzip GCF_000002315.6_GRCg6a_genomic.gff.gz
awk 'BEGIN{OFS="\t";} {if ($3 ~ /exon/) print $1, $4-1, $5}' ${galgal6_gff} | bedtools sort -i - | bedtools merge -i - > galGal6.exon.bed

# get CNEEs
bedtools intersect -v -a galGal6_all_merged.bed -b galGal6.exon.bed > galGal6_final_merged_CNEEs.bed

# add names
awk 'BEGIN {FS="\t"; OFS="\t"} {$4="CNEE"NR; print}' galGal6_final_merged_CNEEs.bed > galGal6_final_merged_CNEEs_named.bed

# replace chr
perl ${replace_chrs_perl} ${trans_matrix} galGal6_final_merged_CNEEs_named.bed > ${final_cnee}
