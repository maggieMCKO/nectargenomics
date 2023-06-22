#!/bin/bash
#SBATCH --time=0:20:00       # it finished in 7 mins
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=4g
#SBATCH --partition=medium
#SBATCH -o splitGffbyChr_%A.out
#SBATCH -e splitGffbyChr_%A.err
#SBATCH -J splitGffbyChr
#SBATCH --mail-type=ALL        # BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mcko@orn.mpg.de


##### SETUP #####
export cneeAna_DIR="/home/mpg08/mko/Nectar/analysis/CNEEanalysis/"
WD="${cneeAna_DIR}00_maf_gff_preprocessing/1_gff/"

export Ref_species="galGal6"
chicken_anno="${cneeAna_DIR}GCF_000002315.6_GRCg6a_genomic.gff" # ori
export split_chicken_anno_DIR="${WD}${Ref_species}_gff_split_galgal/"
mkdir -p ${split_chicken_anno_DIR}

# files needed/generated during fixing chicken anno
trans_matrix=${WD}galgal6a_2315.6_acc3.tsv
fixChr=${chicken_anno%%.gff}_chrFixed2.gff # can be deleted
export fixed=${chicken_anno%%.gff}_chrtabFixed2.gff

##### RUN #####
### 1. FIXING gff chromosome
### 1.1 get corresponding chr matrix (acc) to make trans_matrix
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_assembly_report.txt

sed 's/\r$//g' GCF_000002315.6_GRCg6a_assembly_report.txt | grep -v "^#" | cut -f7,10 | awk '{gsub(/chr/, "chr") ; print}' | awk '{gsub(/\tna/, "\tchrMT") ; print}' > ${trans_matrix}

### 1.2 Fix
perl "${cneeAna_DIR}replace_chrs.pl" ${trans_matrix} ${chicken_anno} > ${fixChr}
awk '{gsub(/Curated\tGenomic/, "Curated_Genomic") ; print}' ${fixChr} > ${fixed}

### 2. split gff by chromosomes [works]
splitGFFbyChr (){
    Chr=$1
    tmp_out="${split_chicken_anno_DIR}${Ref_species}_fixed_${Chr}.gff"
    grep -P "${Chr}\t" ${fixed} > ${tmp_out}
}

export -f splitGFFbyChr

cut -f2 ${trans_matrix} | xargs -n 1 -P ${SLURM_NTASKS} -I {} bash -c 'splitGFFbyChr "$@"' _ {}
