#!/bin/bash

input_ann=/home/mpg08/mko/Nectar/analysis/CNEEanalysis/GCF_000002315.6_GRCg6a_genomic.gff

trans_matrix=galgal6a_2315.6_acc3.tsv
fixChr=${input_ann%%.gff}_chrFixed.gff
fixed=${input_ann%%.gff}_chrtabFixed.gff

# get acc
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_assembly_report.txt
sed 's/\r$//g' GCF_000002315.6_GRCg6a_assembly_report.txt | grep -v "^#" | cut -f7,10 | awk '{gsub(/chr/, "galGal6.chr") ; print}' > ${trans_matrix}

# fix
perl replace_chrs.pl ${trans_matrix} ${input_ann} > ${fixChr}
awk '{gsub(/Curated\tGenomic/, "Curated_Genomic") ; print}' ${fixChr} > ${fixed}

# get maf uniq galgal chr identifiers
# grep "galGal6.chr" multi.anno.maf | awk '{print $2}' | sort -u > maf_galgalchr.txt
