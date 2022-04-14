#!/bin/bash

# WD: "/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/05_perm_phyloacc/1_prepare_new/"
module purge
module load bedtools/2.29.1 ucsc/20160601

cd postPhyloAcc_re2/


# bedtools intersect
intersectFunc (){
    export set=$1
    echo ${set}

    tmp_dir="spatial/"

    # intersect 
    sort -k1,1 -k2,2n "${tmp_dir}${set}_final.cnees.window.sig.bed" | bedtools intersect -a - -b galGal.genes.sorted.bed -loj | sed --expression='s/\.$/0/g' > "${tmp_dir}${set}_final.cnees.window.sig_intersect.bed"

}
export -f intersectFunc

sets=(hummingbirds honeyeaters_pardalote sunbirds_flowerpecker)
parallel -j 4 intersectFunc ::: "${sets[@]}"


intersectFunc2 (){
    export set=$1
    export set2=$2
    echo 'interect2'
    echo -e "set1: ${set}" 
    echo -e "set2: ${set2}" 

    if [ ${set} != ${set2} ]
    then

        echo "do somethin"
        tmp_dir="spatial/"

        sort -k1,1 -k2,2n "${tmp_dir}${set}_final.cnees.window.sig.bed" | bedtools slop -i - -g galGal6.chrom.sizes -b 1000 | bedtools intersect -a - -b "${tmp_dir}${set2}_final.cnees.window.sig.bed" -filenames > "${tmp_dir}${set}_intersect_${set2}.bed"

    fi

}
export -f intersectFunc2

sets=(hummingbirds honeyeaters_pardalote sunbirds_flowerpecker)
parallel -j 4 intersectFunc2 ::: "${sets[@]}" ::: "${sets[@]}"

# awk '$1 == "chr28"' sunbirds_flowerpecker_final.cnees.window.sig_intersect.bed | sed --expression='/RNA/d' | cut -f7 | sort -u

# awk '$8 == "protein_coding"' hummingbirds_final.cnees.window.sig_intersect.bed
# awk '$8 == "protein_coding"' hummingbirds_final.cnees.window.sig_intersect.bed | awk '$1 == "chr23"'  | cut -f7 | sort -u

# awk '$8 == "protein_coding"'  honeyeaters_pardalote_final.cnees.window.sig_intersect.bed | awk '$1 == "chr1"'  | cut -f7 | sort -u