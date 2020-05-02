#!/bin/bash

# trans_matrix=galgal6a_2315.6_acc.tsv
# trans_matrix=galgal6a_2315.6_acc2.tsv
input_ann=GCF_000002315.6_GRCg6a_genomic.gff

fixed=${input_ann%%.gff}_chrtabFixed.gff

# fix
perl replace_chrs.pl ${trans_matrix} ${input_ann} > awk '{gsub(/Curated\tGenomic/, "Curated_Genomic") ; print}' > ${fixed}


# get maf uniq galgal chr identifiers
# grep "galGal6.chr" multi.anno.maf | awk '{print $2}' | sort -u > maf_galgalchr.txt
