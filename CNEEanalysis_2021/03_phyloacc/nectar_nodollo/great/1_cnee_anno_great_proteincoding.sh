#!/bin/bash

# wd: /home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/03_phyloacc/2024ver/nectar_nodollo/great

##### SETUP INPUT and OUTPUT
export cneeAna_DIR="$HOME/Nectar/analysis/CNEEanalysis_2021/"
export chicken_anno="${cneeAna_DIR}00_inputs/GCF_000002315.6_GRCg6a_genomic_chrtabFixed_justchr.gff" # just chr 
export features_bed="${cneeAna_DIR}01b_cnee_prep/bed_outputs/galGal6_final_conserved_CNEEs.bed"


# I. createRegulatoryDomains
# 1. make chrom.sizes
## based on script: 1_cnee_shuf_spaital.sh
export ucsc_path=$HOME/Tools/ucsc_liftover/bin/
${ucsc_path}faToTwoBit "${cneeAna_DIR}/00_inputs/genomes/galGal6.fa" galGal6.fa.2bit
${ucsc_path}twoBitInfo galGal6.fa.2bit stdout | sort -k2rn > galGal6.chrom.sizes


# 2. create TSS.in [4 cols: chromosome, transcription_start_site, strand, geneName]
gene_bed="galGal6_gene_tss.bed"
awk '$3 == "gene"'  ${chicken_anno} | awk -F "\t" '{OFS=FS}{split($9, a, ";"); $9=a[1]; print $1,$2,$3,$4,$5,$6,$7,$8,$9 }' | awk -F "\t" '{OFS=FS}{split($9, a, "gene-"); $9=a[2]; print $1,$2,$3,$4,$5,$6,$7,$8,$9 }' | awk '{ if ($1 != "chrMT") { print } }' | gff2bed | awk 'BEGIN { OFS = "\t" }  {print $1,$2,$6,$10}' | sort -k1,1V -k2,2n > ${gene_bed} 

# way3
grep protein_coding ${chicken_anno} > galGal.genes.bed # 17485

module load rev/23.12 anaconda3/2023.09-0
source activate R432 

export wd="/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/03_phyloacc/2024ver/nectar_nodollo/great/"
Rscript 1b_proteincoding_subset.R ${wd}${gene_bed} ${wd}galGal.genes.bed "${wd}galGal6_gene_tss_proteincoding.bed"

# 3. run createRegulatoryDomains
	# $ createRegulatoryDomains  -maxExtension -basalUpstream -basalDownstream TSS.in chrom.sizes oneClosest regDoms.out
	# 		1. TSS.in
	# 		2. chrom.sizes
	# 		3. oneClosest|twoClosest|basalPlusExtension
	# 		4. regDoms.out

module purge
module load rev/11.06 GREAT/1.5 

gene_bed="galGal6_gene_tss_proteincoding.bed"
createRegulatoryDomains ${gene_bed} galGal6.chrom.sizes basalPlusExtension galGal6_gene_tss_great_basalPlusExtension_proteincoding.bed
# wc -l galGal6_gene_tss_great_basalPlusExtension.bed # 24108
# wc -l galGal6_gene_tss_great_basalPlusExtension_proteincoding.bed # 17472
# chromStart and chromEnd of the 4 outputs (a bed file) are different


# II. create antigap.bed for running calculateBinomialP
# https://bioinformatics.stackexchange.com/questions/3526/bed-file-with-n-regions-of-grch38-reference
# twoBitInfo file.fa -nBed output.bed

# bed for Ns
${ucsc_path}twoBitInfo galGal6.fa.2bit -nBed galGal6_gap.bed 

# gff to bed [ used the same gff as those were used to create regulatory domain ]
module purge
module load rev/20.12 bedtools/2.29.1
# ls /usr/product/bioinfo/SL_7.0/BIOINFORMATICS/UCSC/20160601

awk -F "\t" '{OFS=FS}{split($9, a, ";"); $9=a[1]; print $1,$2,$3,$4,$5,$6,$7,$8,$9 }' ${chicken_anno} | awk -F "\t" '{OFS=FS}{split($9, a, "gene-"); $9=a[2]; print $1,$2,$3,$4,$5,$6,$7,$8,$9 }' | awk '{ if ($1 != "chrMT") { print } }' | gff2bed | awk 'BEGIN { OFS = "\t" }  {print $1,$2,$3,$4,$5,$6,$7,$8,$9}' | sort -k1,1V -k2,2n | bedtools merge > galGal6_gff2bed.bed

# genome - substract gap
bedtools subtract -a galGal6_gff2bed.bed -b galGal6_gap.bed | awk '{print $0, "."}' > galGal6_antigap.bed 