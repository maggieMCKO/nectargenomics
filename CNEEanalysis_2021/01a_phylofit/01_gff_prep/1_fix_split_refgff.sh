#!/bin/bash
#SBATCH --time=0:10:00       # it finished in 7 mins
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=4g
#SBATCH --partition=medium
#SBATCH -o splitChr_%A.out
#SBATCH -e splitChr_%A.err
#SBATCH -J splitChr
#SBATCH --mail-type=ALL        # BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mcko@orn.mpg.de

# srun -p int --pty -t 2:00:00 -n 1 bash

##### SETUP #####

export cneeAna_DIR="$HOME/Nectar/analysis/CNEEanalysis_2021/"
export WD="${cneeAna_DIR}01a_phylofit/01_gff_prep/"

export Ref_species="galGal6"
export chicken_anno="${cneeAna_DIR}00_inputs/GCF_000002315.6_GRCg6a_genomic.gff" # ori

export split_chicken_anno_DIR="${WD}${Ref_species}_gff_split/"
mkdir -p ${split_chicken_anno_DIR}

# files needed/generated during fixing chicken anno
fixChr=${chicken_anno%%.gff}_chrFixed2.gff # can be deleted
export fixed=${chicken_anno%%.gff}_chrtabFixed2.gff

fixChr_justchr=${chicken_anno%%.gff}_chrFixed_justchr.gff # can be deleted
export fixed_justchr=${chicken_anno%%.gff}_chrtabFixed_justchr.gff


##### RUN #####
### 1. FIXING gff chromosome
# 1.1 download galgal assembly report
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_assembly_report.txt

## make acc and chromosome converion table (galGal6.chr): for spliting gff
trans_matrix_for_gff=galgal6a_2315.6_acc_galgal6chr_matrix.tsv
sed 's/\r$//g' GCF_000002315.6_GRCg6a_assembly_report.txt | grep -v "^#" | cut -f7,10 | awk '{gsub(/chr/, "galGal6.chr") ; print}' | awk '{gsub(/\tna/, "\tchrMT") ; print}' > ${trans_matrix_for_gff}

## make acc and chromosome converion table (.chr): for spliting maf
trans_matrix_for_maf=galgal6a_2315.6_acc_chr_matrix.tsv
sed 's/\r$//g' GCF_000002315.6_GRCg6a_assembly_report.txt | grep -v "^#" | cut -f7,10 | awk '{gsub(/chr/, "chr") ; print}' | awk '{gsub(/\tna/, "\tchrMT") ; print}' > ${trans_matrix_for_maf}

### 1.2 Fix
perl replace_chrs.pl ${trans_matrix_for_gff} ${chicken_anno} > ${fixChr}
awk '{gsub(/Curated\tGenomic/, "Curated_Genomic") ; print}' ${fixChr} > ${fixed}

perl replace_chrs.pl ${trans_matrix_for_maf} ${chicken_anno} > ${fixChr_justchr}
awk '{gsub(/Curated\tGenomic/, "Curated_Genomic") ; print}' ${fixChr_justchr} > ${fixed_justchr}

### 2. split gff by chromosomes
splitGFFbyChr (){
    Chr=$1
    tmp_out="${split_chicken_anno_DIR}${Ref_species}_fixed_${Chr}.gff"
    grep -P "${Chr}\t" ${fixed} > ${tmp_out}
}

export -f splitGFFbyChr

cut -f2 ${trans_matrix_for_gff} | parallel splitGFFbyChr {}
