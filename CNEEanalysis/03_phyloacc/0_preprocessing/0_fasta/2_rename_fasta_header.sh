#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --ntasks=30
#SBATCH --mem=60g
#SBATCH --partition=medium
#SBATCH -o faRename_%A.out
#SBATCH -e faRename_%A.err
#SBATCH -J faRename
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mcko@orn.mpg.de
#SBATCH -C "scratch|scratch2"

export cneeAna_DIR="/home/mpg08/mko/Nectar/analysis/CNEEanalysis/"

export WD="${cneeAna_DIR}03_phyloacc/preprocessing/0_fasta/"
export specieslist="${cneeAna_DIR}03_phyloacc/preprocessing/1_split_sp/Sp.list"

export maf_DIR="${cneeAna_DIR}00_maf_gff_preprocessing/2_maf/"
export input_maf="${maf_DIR}multi.anno_add1st2.maf"

export genome_fa_DIR="${cneeAna_DIR}03_phyloacc/preprocessing/0_fasta/ori/"
export assembly_report_DIR="${cneeAna_DIR}03_phyloacc/preprocessing/0_fasta/assembly_report/"

####### RUN
### create species list
mafSpeciesList in.maf out.lst
mafSpeciesList ${input_maf} ${specieslist}


splitsp (){
    SP=$1

    # Setup rename fa file names
    tmp_fa="${genome_fa_DIR}${SP}.fasta" # from
    Renamed_DIR="${cneeAna_DIR}03_phyloacc/preprocessing/0_fasta/renamed_fa/"
    mkdir -p ${Renamed_DIR}
    tmp_fa_renamed="${Renamed_DIR}${SP}_renamed.fasta"

    if [ "$SP" = "HLacaPus1" ] || [ "$SP" = "HLfloFus1" ] || [ "$SP" = "HLfurRuf1" ] || [ "$SP" = "HLphyNov1" ] || [ "$SP" = "HLtriMol2" ] ; then

        echo "current SP: $SP - OURS"

        # rename serial numbers
        awk -v species="$SP" '/^>/{printf ">%s_%.8d\n", species, ++i ; next}{print}'  ${tmp_fa} > ${tmp_fa_renamed}

    elif [ "$SP" = "HLcolLiv2" ] || [ "$SP" = "HLfalTin1" ] || [ "$SP" = "HLlicCas1" ] || [ "$SP" = "HLmalCya1" ] || [ "$SP" = "HLzosLat1" ] || [ "$SP" = "HLchaPel1" ] || [ "$SP" = "HLtaeGut4" ] ; then
        echo "current SP: $SP - NCBI"
        echo "current SP: $SP - shorten header"

        # rename
        awk '{
            if ($1 ~ /^>/) {
                split($0,a," "); split(a[1], b, "\."); print b[1]
            } else {
                print $0 }
        }' ${tmp_fa} > ${tmp_fa_renamed}


    else
        echo "current SP: $SP - NCBI"

        if [ "$SP" = "galGal6" ]; then
            echo "current SP: $SP - col 7 to col 10"

            # make matrix
            tmp_assembly_report="${assembly_report_DIR}${SP}_assembly_report.txt"
            tmp_matrix="${assembly_report_DIR}${SP}_matrix.txt"
            sed 's/\r$//g' ${tmp_assembly_report} | grep -v "^#" | cut -f7,10 | awk '{gsub(/\tna/, "\tchrMT") ; print}' > ${tmp_matrix}
            # | awk '{gsub(/chr/, "chr") ; print}' |

        elif [ "$SP" = "ficAlb2" ] ; then
            echo "current SP: $SP - col 7 to col 1"

            # make matrix
            tmp_assembly_report="${assembly_report_DIR}${SP}_assembly_report.txt"
            tmp_matrix="${assembly_report_DIR}${SP}_matrix.txt"

            # cut -f 15  ../1_way2/CneeMappedToGg6_PSL_DIR/ficAlb2.psl | sort -u

            sed 's/\r$//g' ${tmp_assembly_report} | grep -v "^#" |
            awk -F "\t" '{OFS=FS}{split($5, a, "\."); $5=a[1]; print }' |
            awk -F "\t" '{OFS=FS}{
                if ($1 ~ /^N/) {$5=$5"_random"; print } else {print $0}
            }' |
            awk -F "\t" '{OFS=FS}{
                if ($1 ~ /chr/) {print $7, $1}
            else { if ($1 ~ /LG/) {print $7, "chr"$3 } else {print $7, "chr"$3"_"$5 } }
            }' |
            awk -F "\t" '{OFS=FS}{gsub(/chrna/, "chrUn") ; print}' |
            awk -F "\t" '{OFS=FS}{gsub(/\tna/, "\tchrMT") ; print}'  > ${tmp_matrix}

        else
            echo "current SP: $SP - col 7 to col 5"
            # make matrix
            tmp_assembly_report="${assembly_report_DIR}${SP}_assembly_report.txt"
            tmp_matrix="${assembly_report_DIR}${SP}_matrix.txt"
            sed 's/\r$//g' ${tmp_assembly_report} | grep -v "^#" | awk -F "\t" '{OFS=FS}{print $7, $5}' | awk -F "\t" '{OFS=FS}{gsub(/\tna/, "\tchrMT") ; print}' | awk -F "\t" '{OFS=FS}{split($2, a, "\."); print $1, a[1]}' > ${tmp_matrix}

        fi

        # rename
        awk 'FNR==NR {
            hash[">"$1]=$2;next
        }
        {
        if ($1 ~ /^>/)
            print ">"hash[$1];
        else
            print $1;
        }' ${tmp_matrix} ${tmp_fa} > ${tmp_fa_renamed}
    fi


}

export -f splitsp

cut -f1 ${specieslist} | xargs -n 1 -P ${SLURM_NTASKS} -I {} bash -c 'splitsp "$@"' _ {}
