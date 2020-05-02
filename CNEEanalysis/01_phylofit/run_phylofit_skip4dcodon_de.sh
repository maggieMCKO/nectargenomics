#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=4g
#SBATCH --partition=medium
#SBATCH -o phylofit_%A.out
#SBATCH -e phylofit_%A.err
#SBATCH -J phylofit
#SBATCH --mail-type=ALL        # BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mcko@orn.mpg.de

##### SETUP INPUT and OUTPUT
WD="/home/mpg08/mko/Nectar/analysis/CNEEanalysis/"
input_maf="${WD}multi.anno.maf"
input_tree="${WD}NectarTree01.nw"
chicken_anno="${WD}GCF_000002315.6_GRCg6a_genomic_chrtabFixed_CDS.gff" # chr need to match maf's chr

out_DIR="${WD}01_phylofit/01_msa/"
mkdir -p ${out_DIR}
codons_4d="${out_DIR}4d-codons.ss"
sites_4d="${out_DIR}4d-sites.ss"

out_DIR_4d="${WD}nonconserved_4d"


##### SETUP programs
phastPATH="/home/mpg08/mko/Tools/phast/bin/"
phyloFit_PATH="${phastPATH}phyloFit"
msa_view_PATH="${phastPATH}msa_view"


##### RUN
# extract 4d codons from an alignment
# ${msa_view_PATH} ${input_maf} --4d --features ${chicken_anno} > ${codons_4d}

# extract 4d sites from an alignment
${msa_view_PATH} ${input_maf} --in-format MAF --out-format SS --tuple-size 1 > ${sites_4d}

# estimate A nonconserved phylogenetic model
${phyloFit_PATH} --tree ${input_tree} --msa-format SS --out-root ${out_DIR_4d} ${sites_4d}
