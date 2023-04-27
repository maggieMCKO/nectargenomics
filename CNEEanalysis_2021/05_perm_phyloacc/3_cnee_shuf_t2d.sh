#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --ntasks=9
#SBATCH --mem-per-cpu=5g
#SBATCH --partition=fat
#SBATCH -o log1_t2dperm_go_%A.o
#SBATCH -e log1_t2dperm_go_%A.e
#SBATCH -J perm_go


# modified from https://github.com/tsackton/ratite-genomics/blob/master/04_wga/02_ce_id/postprocess_ces.sh


# WD: "/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/05_perm_phyloacc/1_rmPar/"
# 20230427


source 0_setup_paths.sh
module purge
module load bedtools/2.29.1 anaconda2/2019.10
# source activate R_phylo


mkdir -p postPhyloAcc
cd postPhyloAcc/

### 3. input generation for GO permutations 
# nb: could also write quick script (like replace_chr.pl) to replace the gene symbols with NCBI gene IDs in cnee_perms.counts.bed

mkdir -p t2d_perms

# # pull out gene symbols & NCBI IDs
# cat galGal.gff.bed | python3 ../gene_ncbi_names.py | awk '{ if ($4 != "") { print } }' > galGal.ncbigenes.bed

# # sort 
# bedtools sort -i galGal.ncbigenes.bed > galGal.ncbigenes.sorted.bed

# # set 100 kb window around genes
# bedtools slop -i galGal.ncbigenes.sorted.bed -g galGal6.chrom.sizes_try2 -b 100000 > galGal.ncbislop.bed

bedtools slop -i ../ggal6/galGal6.ncbigenes.sorted.bed -g galGal6.chrom.sizes_try2 -b 100000 > galGal.ncbislop.bed # new: ../ggal6/galGal6.ncbigenes.sorted.bed

# annotate slopped BED with accel CNEEs 
export output_dir=/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/04_post_phyloacc/outputs_rmPar/

export sets=(4sets honeyeaters hummingbirds parrots sunbirds_flowerpecker)
# export t2d_sets=(t2d_strong t2d_intermediate t2d_all )

for set in "${sets[@]}"; do
	echo $set
	# t2d_strong
	bedtools annotate -i galGal.ncbislop.bed -files "${output_dir}Target_NeFr40_thre_0.9_logBF3_10_${set}.bed" "${output_dir}Target_NeFr40_thre_0.9_logBF3_10_t2d_strong.bed" galGal6_final_conserved_CNEEs_sorted.bed -counts > "t2d_perms/cnee_ncbigene100kb_${set}_t2d_strong_acc.bed" 
	awk '$8 > 0 {print}' "t2d_perms/cnee_ncbigene100kb_${set}_t2d_strong_acc.bed" | sed 's/ /\t/g' - > "t2d_perms/cnee_ncbigene100kb_${set}_t2d_strong_acc_clean.bed"


	# t2d_intermediate
	bedtools annotate -i galGal.ncbislop.bed -files "${output_dir}Target_NeFr40_thre_0.9_logBF3_10_${set}.bed" "${output_dir}Target_NeFr40_thre_0.9_logBF3_10_t2d_intermediate.bed" galGal6_final_conserved_CNEEs_sorted.bed -counts > "t2d_perms/cnee_ncbigene100kb_${set}_t2d_intermediate_acc.bed" 
	awk '$8 > 0 {print}' "t2d_perms/cnee_ncbigene100kb_${set}_t2d_intermediate_acc.bed" | sed 's/ /\t/g' - > "t2d_perms/cnee_ncbigene100kb_${set}_t2d_intermediate_acc_clean.bed"


	# t2d_all (weak)
	bedtools annotate -i galGal.ncbislop.bed -files "${output_dir}Target_NeFr40_thre_0.9_logBF3_10_${set}.bed" "${output_dir}Target_NeFr40_thre_0.9_logBF3_10_t2d_all.bed" galGal6_final_conserved_CNEEs_sorted.bed -counts > "t2d_perms/cnee_ncbigene100kb_${set}_t2d_all_acc.bed" 
	awk '$8 > 0 {print}' "t2d_perms/cnee_ncbigene100kb_${set}_t2d_all_acc.bed" | sed 's/ /\t/g' - > "t2d_perms/cnee_ncbigene100kb_${set}_t2d_all_acc_clean.bed"

done


# do permutations
go_perms_func (){
	export set=$1
	export n=$2

	echo ${set}
	echo ${n}

	# do permutations
	tmp_dir="t2d_perms/perms_${set}/"
	mkdir -p ${tmp_dir}
	for i in {0001..1000}; 
	do
	  shuf -n ${n} galGal6_final_conserved_CNEEs_sorted.bed > "${tmp_dir}shuffle${i}.bed"
	done
	# annotate 

	# t2d_strong
	bedtools annotate -i "t2d_perms/cnee_ncbigene100kb_${set}_t2d_strong_acc_clean.bed" -files ${tmp_dir}*.bed -counts > "t2d_perms/cnee_goperms_${set}_t2d_strong.counts.bed"

	# t2d_intermediate
	bedtools annotate -i "t2d_perms/cnee_ncbigene100kb_${set}_t2d_intermediate_acc_clean.bed" -files ${tmp_dir}*.bed -counts > "t2d_perms/cnee_goperms_${set}_t2d_intermediate.counts.bed"

	# t2d_all (weak)
	bedtools annotate -i "t2d_perms/cnee_ncbigene100kb_${set}_t2d_all_acc_clean.bed" -files ${tmp_dir}*.bed -counts > "t2d_perms/cnee_goperms_${set}_t2d_all.counts.bed"
}
export -f go_perms_func
n_perms=(3269 1691 1125 481) # num of cnees
# n_perms=(159 589 1408) # cant use this script

parallel --xapply -j 4 go_perms_func ::: "${sets[@]}" ::: "${n_perms[@]}"
