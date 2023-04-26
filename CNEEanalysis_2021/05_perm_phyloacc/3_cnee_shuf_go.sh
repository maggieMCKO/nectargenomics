#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --ntasks=6
#SBATCH --mem-per-cpu=5g
#SBATCH --partition=fat
#SBATCH -o log1_perm_go_%A.o
#SBATCH -e log1_perm_go_%A.e
#SBATCH -J perm_go


# modified from https://github.com/tsackton/ratite-genomics/blob/master/04_wga/02_ce_id/postprocess_ces.sh


# WD: "/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/05_perm_phyloacc/1_rmPar/"

source 0_setup_paths.sh
module purge
module load bedtools/2.29.1 anaconda2/2019.10
# source activate R_phylo


mkdir -p postPhyloAcc
cd postPhyloAcc/

### 3. input generation for GO permutations 
# nb: could also write quick script (like replace_chr.pl) to replace the gene symbols with NCBI gene IDs in cnee_perms.counts.bed

mkdir -p go_perms

# # pull out gene symbols & NCBI IDs
# cat galGal.gff.bed | python3 ../gene_ncbi_names.py | awk '{ if ($4 != "") { print } }' > galGal.ncbigenes.bed

# # sort 
# bedtools sort -i galGal.ncbigenes.bed > galGal.ncbigenes.sorted.bed

# # set 100 kb window around genes
# bedtools slop -i galGal.ncbigenes.sorted.bed -g galGal6.chrom.sizes_try2 -b 100000 > galGal.ncbislop.bed

bedtools slop -i ../ggal6/galGal6.ncbigenes.sorted.bed -g galGal6.chrom.sizes_try2 -b 100000 > galGal.ncbislop.bed # new: ../ggal6/galGal6.ncbigenes.sorted.bed

# annotate slopped BED with accel CNEEs 
export output_dir=/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/04_post_phyloacc/outputs_rmPar/

export sets=(4sets 2wayconvergent 3wayconvergent honeyeaters hummingbirds parrots sunbirds_flowerpecker)
# export sets=(genelv_conv_4way genelv_conv_3way genelv_conv_2way ) # cant use this script

for set in "${sets[@]}"; do
	echo $set
	bedtools annotate -i galGal.ncbislop.bed -files "${output_dir}Target_NeFr40_thre_0.9_logBF3_10_${set}.bed" galGal6_final_conserved_CNEEs_sorted.bed -counts > "go_perms/cnee_ncbigene100kb_${set}_acc.bed" 
	awk '$8 > 0 {print}' "go_perms/cnee_ncbigene100kb_${set}_acc.bed" | sed 's/ /\t/g' - > "go_perms/cnee_ncbigene100kb_${set}_acc_clean.bed"
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
	bedtools annotate -i "go_perms/cnee_ncbigene100kb_${set}_acc_clean.bed" -files ${tmp_dir}*.bed -counts > "go_perms/cnee_goperms_${set}.counts.bed"
}
export -f go_perms_func
n_perms=(3269 477 37 523 1691 1125 481) # num of cnees
# n_perms=(159 589 1408) # cant use this script

parallel --xapply -j 4 go_perms_func ::: "${sets[@]}" ::: "${n_perms[@]}"
