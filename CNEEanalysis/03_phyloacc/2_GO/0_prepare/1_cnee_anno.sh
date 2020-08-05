#!/bin/bash

# modified from https://github.com/tsackton/ratite-genomics/blob/master/04_wga/02_ce_id/postprocess_ces.sh

source setup_paths.sh

module purge
module load BEDTOOLS UCSC/20160601
# ls /usr/product/bioinfo/SL_7.0/BIOINFORMATICS/UCSC/20160601

export WD="/home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/2_GO/0_prepare/"

# 1. get galGal6_gene.bed, galGal6_exon.bed, galGal6_cds.bed
tmp_filename=$(basename ${chicken_anno})

# gene
gene_gff="${WD}galGal6_gene.gff"
gene_bed="${WD}galGal6_gene.bed"
awk '{if ($3 ~ /gene/) print $0}' ${chicken_anno} | awk -F "\t" '{OFS=FS}{split($9, a, ";"); $9=a[1]; print }' | awk -F "\t" '{OFS=FS}{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > ${gene_gff}
gff2bed < ${gene_gff} > ${gene_bed} # from bedops pacakge


# exon [seems not necessary]
exon_gff="${WD}galGal6_exon.gff"
exon_bed="${WD}galGal6_exon.bed"
awk '{if ($3 ~ /exon/) print $0}' ${chicken_anno} | awk -F "\t" '{OFS=FS}{split($9, a, ";"); $9=a[1]; print }' | awk -F "\t" '{OFS=FS}{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > ${exon_gff}
bedtools sort -i ${exon_gff} | bedtools merge -i - -s -d -1 -c 4 -o distinct  > ${exon_bed}

# cds [seems not necessary]
cds_gff="${WD}galGal6_cds.gff"
cds_bed="${WD}galGal6_cds.bed"
awk '{if ($3 ~ /CDS/) print $0}' ${chicken_anno} | awk -F "\t" '{OFS=FS}{split($9, a, ";"); $9=a[1]; print }' | awk -F "\t" '{OFS=FS}{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > ${cds_gff}
bedtools sort -i ${cds_gff} | bedtools merge -i - -s -d -1 -c 4 -o distinct  > ${cds_bed}

# sort CNEE, otherwise, get error when using bedtools closest
features_bed_sorted="${WD}galGal6_final_merged_CNEEs_named_fixchr_justchr_sorted.bed"
sort -k1,1 -k2,2n ${features_bed} > ${features_bed_sorted}

# A get closest gene
cloest_gene="${WD}galGal6_final_merged_CNEEs_cloest_genes.bed"
bedtools closest -a ${features_bed_sorted} -b ${gene_bed} -D "b" -t "all" > ${cloest_gene}

# B get closest genes defined by approximate TSS, sort, otherwise, get error when using bedtools closest
gene_TSS_bed="${WD}galGal6_approxTSS.bed"
./get_approx_TSS.sh ${gene_bed} | sort -k1,1 -k2,2n > ${gene_TSS_bed}

cloest_TSS="${WD}galGal6_final_merged_CNEEs_cloest_TSS.bed"
bedtools closest -a ${features_bed_sorted} -b ${gene_TSS_bed} -D "b" -t "all" > ${cloest_TSS}
