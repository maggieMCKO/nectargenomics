#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=4g
#SBATCH --partition=medium
#SBATCH -o phylofit_%A.out
#SBATCH -e phylofit_%A.err
#SBATCH -J phylofit
#SBATCH --mail-type=ALL        # BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mcko@orn.mpg.de

##### SETUP INPUT and OUTPUT
export cneeAna_DIR="/home/mpg08/mko/Nectar/analysis/CNEEanalysis/"
export input_tree="${cneeAna_DIR}NectarTree01.nw"
export chicken_anno="${cneeAna_DIR}GCF_000002315.6_GRCg6a_genomic_chrtabFixed2.gff" # chr need to match maf's chr with species (galgal.chrx)

export maf_DIR="/home/mpg08/mko/Nectar/analysis/CNEEanalysis/00_maf_gff_preprocessing/2_maf/"
export input_maf="${maf_DIR}multi.anno.maf"
export maf_processed_DIR="${maf_DIR}Processed"

export gff_DIR="${cneeAna_DIR}00_maf_gff_preprocessing/1_gff/"
export trans_matrix="${gff_DIR}galgal6a_2315.6_acc3.tsv"
export gff_processed_DIR="${gff_DIR}galGal6_gff_split_galgal/"

export WD="${cneeAna_DIR}01_phylofit/"
export out_DIR="${WD}01_msa3/"
mkdir -p ${out_DIR}

sp_list="galGal6,melUnd1,HLtriMol2,HLstrHab1,falPer1,falChe1,HLfalTin1,HLtytAlb2,halLeu1,aptFor1,HLcolLiv2,HLempTra1,HLcalPug1,cucCan1,HLtaeGut4,HLlicCas1,HLphyNov1,HLacaPus1,HLfurRuf1,ficAlb2,HLparMaj1,HLserCan1,HLmalCya1,pseHum1,HLzosLat1,HLcorCor3,HLcalAnn5,HLfloFus1,HLchaPel1,opiHoa1,HLphaSup1"
concat_4ds="${out_DIR}conca_4dSites.ss"
out_DIR_4d="${out_DIR}nonconserved_4d"

##### SETUP programs
export phastPATH="/home/mpg08/mko/Tools/phast/bin/"
export phyloFit_PATH="${phastPATH}phyloFit"
export msa_view_PATH="${phastPATH}msa_view"


##### RUN
### 1. extract 4d codons and 4d sites from an alignment (by chromosome)
Get4dCodonSitebyChr (){
    Chr=$1
    chr_clean=$(echo $Chr | sed 's/galGal6.//')
    # echo $Chr

    # input maf
    chr_maf="${chr_clean}.maf"
    tmp_maf=$(find ${maf_processed_DIR} | grep ${chr_maf})
    line_maf=$(wc -l $tmp_maf)
    echo -e "maf: ${line_maf}"
    line_maf_cut=$(echo ${line_maf} | cut -d ' ' -f1 )

    # input gff
    chr_gff="galGal6_fixed_${Chr}.gff"
    tmp_gff=$(find ${gff_processed_DIR} | grep ${chr_gff})
    line_gff=$(wc -l $tmp_gff)
    echo -e "gff: ${line_gff}"
    line_gff_cut=$(echo ${line_gff} | cut -d ' ' -f1 )

    if (( $line_maf_cut > 2 )) && (( $line_gff_cut > 1)); then
        echo $chr_clean

        # output 4d codons and sites
        codons_4d="${out_DIR}${chr_clean}_4dCodons.ss"
        sites_4d="${out_DIR}${chr_clean}_4dSites.ss"

        # 4d codon
        ${msa_view_PATH} ${tmp_maf} --4d --features ${tmp_gff} > ${codons_4d}

        # 4d sites
        ${msa_view_PATH} ${codons_4d} --in-format SS --out-format SS --tuple-size 1 > ${sites_4d}
    fi
}

export -f Get4dCodonSitebyChr
cut -f2 ${trans_matrix} | xargs -n 1 -P 1 -I {} bash -c 'Get4dCodonSitebyChr "$@"' _ {}

### 2. concatenate 4d sites
${msa_view_PATH} --aggregate ${sp_list} --in-format SS --out-format SS --unordered-ss chr*_4dSites.ss > ${concat_4ds}

### 3. run phyloFit: estimate A nonconserved phylogenetic model
${phyloFit_PATH} --tree ${input_tree} --msa-format SS --out-root ${out_DIR_4d} ${concat_4ds}
