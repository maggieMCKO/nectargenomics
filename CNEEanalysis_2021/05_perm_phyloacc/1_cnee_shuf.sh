#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --ntasks=24
#SBATCH --mem-per-cpu=5g
#SBATCH --partition=medium
#SBATCH -o log1_perm_%A.o
#SBATCH -e log1_perm_%A.e
#SBATCH -J perm


# modified from https://github.com/tsackton/ratite-genomics/blob/master/04_wga/02_ce_id/postprocess_ces.sh


# WD: "/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/05_perm_phyloacc/1_prepare/"

source 0_setup_paths.sh
module purge
module load bedtools/2.29.1 anaconda2/2019.10
source activate R_phylo


mkdir -p postPhyloAcc_re2
cd postPhyloAcc_re2/

### 1. spatial enrichment analyses input generation

mkdir -p spatial


# sort full final list
bedtools sort -i ${features_bed} > galGal6_final_conserved_CNEEs_sorted.bed

# make chrom.sizes
export ucsc_path=$HOME/Tools/ucsc_liftover/bin/
${ucsc_path}faToTwoBit "${cneeAna_DIR}/00_inputs/genomes/galGal6.fa" stdout | ${ucsc_path}twoBitInfo stdin stdout | sort -k2,2nr | awk '{ if ($2 > 0) { print } }'  > galGal6.chrom.sizes

# create 100kb windows with 50kb slide
bedtools makewindows -g galGal6.chrom.sizes -w 100000 -s 50000 > galGal6.windows.bed

# bin total CNEEs list into windows
bedtools intersect -a galGal6.windows.bed -b galGal6_final_conserved_CNEEs_sorted.bed -loj | cut -f1,2,3,7 | sed --expression='s/\.$/0/g' > window.cnees.bed

# use output from 04_PhyloAcc/07_phyloP_cleanup_cnees.R to bin accelerated CNEEs list into windows
output_dir=/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/04_post_phyloacc/outputs/
sets=(hummingbirds parrots honeyeaters_pardalote sunbirds_flowerpecker)
for set in "${sets[@]}"; do
	echo $set
	bedtools intersect -a galGal6.windows.bed -b "${output_dir}Target_NeFr40_thre_0.9_logBF3_10_${set}.bed" -loj | cut -f1,2,3,7 | sed --expression='s/\.$/0/g' > spatial/window.acc.cnees_${set}.bed
done


### 2. input generation for the assessment of genes with evidence for excess of nearby accelerated CNEEs

mkdir -p gene


# gff to bed
# awk '{print NF}' ${chicken_anno} | sort -u 
awk '$3 == "gene"' ${chicken_anno} | awk '{gsub(/chrMT/, "chrM") ; print}' | awk '$1 != "chrM"' | sed --expression='/gene_biotype=tRNA/d' | sed --expression='/gene_biotype=rRNA/d' > galGal.genes.gff
gff2bed < galGal.genes.gff | bedtools sort -i - > galGal.gff.bed # gff2bed: from bedops pacakge 
# no dup
# awk '{print NF}' galGal.gff.bed | sort -u 


# pull out genes
cat galGal.gff.bed | python3 ../genenames.py  > galGal.genes.bed
# dup happened here
# awk '{print NF}' galGal.genes.bed | sort -u 

# sort 
bedtools sort -i galGal.genes.bed > galGal.genes.sorted.bed

# set 100 kb window around genes
bedtools slop -i galGal.genes.sorted.bed -g galGal6.chrom.sizes -b 100000 > galGal.slop.bed 

# annotate slopped BED with accel CNEEs 
# output_dir=/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/04_post_phyloacc/outputs/
# sets=(hummingbirds parrots honeyeaters_pardalote sunbirds_flowerpecker)
for set in "${sets[@]}"; do
	echo $set
	bedtools annotate -i galGal.slop.bed -files "${output_dir}Target_NeFr40_thre_0.9_logBF3_10_${set}.bed" galGal6_final_conserved_CNEEs_sorted.bed -counts > "gene/cnee_gene100kb_${set}_acc.bed" 
	awk '$7 > 0 {print}' "gene/cnee_gene100kb_${set}_acc.bed" | sed 's/ /\t/g' - > "gene/cnee_gene100kb_${set}_acc_clean.bed"
done

# do permutations
cnee_perms_func (){
	export set=$1
	export n=$2

	echo ${set}
	echo ${n}

	# do permutations
	tmp_dir="gene/cnee_perms_${set}/"
	mkdir -p ${tmp_dir}
	for i in {0001..1000}; 
	do
	  shuf -n ${n} galGal6_final_conserved_CNEEs_sorted.bed > "${tmp_dir}shuffle${i}.bed"
	done
	# annotate 
	bedtools annotate -i "gene/cnee_gene100kb_${set}_acc_clean.bed" -files ${tmp_dir}*.bed -counts > "gene/cnee_perms_${set}.counts.bed"
}
export -f cnee_perms_func
sets=(hummingbirds parrots honeyeaters_pardalote sunbirds_flowerpecker)
n_perms=(1691 1125 655 481)

parallel --xapply -j 4 cnee_perms_func ::: "${sets[@]}" ::: "${n_perms[@]}"


### 3. input generation for GO permutations 
# nb: could also write quick script (like replace_chr.pl) to replace the gene symbols with NCBI gene IDs in cnee_perms.counts.bed

mkdir -p go_perms

# pull out gene symbols & NCBI IDs
cat galGal.gff.bed | python3 ../gene_ncbi_names.py | awk '{ if ($4 != "") { print } }' > galGal.ncbigenes.bed

# sort 
bedtools sort -i galGal.ncbigenes.bed > galGal.ncbigenes.sorted.bed

# set 100 kb window around genes
bedtools slop -i galGal.ncbigenes.sorted.bed -g galGal6.chrom.sizes -b 100000 > galGal.ncbislop.bed

# annotate slopped BED with accel CNEEs 

for set in "${sets[@]}"; do
	echo $set
	bedtools annotate -i galGal.ncbislop.bed -files "${output_dir}Target_NeFr40_thre_0.9_logBF3_10_${set}.bed" galGal6_final_conserved_CNEEs_sorted.bed -counts > "go_perms/cnee_ncbigene100kb_${set}_acc.bed" 
	awk '$9 > 0 {print}' "go_perms/cnee_ncbigene100kb_${set}_acc.bed" | sed 's/ /\t/g' - > "go_perms/cnee_ncbigene100kb_${set}_acc_clean.bed"
done


# do permutations
go_perms_func (){
	export set=$1
	export n=$2

	echo ${set}
	echo ${n}

	# do permutations
	tmp_dir="go_perms/perms_${set}/"
	mkdir -p ${tmp_dir}
	for i in {0001..1000}; 
	do
	  shuf -n ${n} galGal6_final_conserved_CNEEs_sorted.bed > "${tmp_dir}shuffle${i}.bed"
	done
	# annotate 
	bedtools annotate -i "go_perms/cnee_ncbigene100kb_${set}_acc.bed" -files ${tmp_dir}*.bed -counts > "go_perms/cnee_goperms_${set}.counts.bed"
}
export -f go_perms_func
sets=(hummingbirds parrots honeyeaters_pardalote sunbirds_flowerpecker)
n_perms=(1691 1125 655 481)

parallel --xapply -j 4 go_perms_func ::: "${sets[@]}" ::: "${n_perms[@]}"
