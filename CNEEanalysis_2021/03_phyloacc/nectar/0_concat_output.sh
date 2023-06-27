# https://github.com/sjswuitchik/duck_comp_gen/blob/master/03a_cnee_analysis/06_PhyloAcc/concat_output.sh
# concat_output.sh
# used in parse_all.sh

awk '
    FNR==1 && NR!=1 { while (/^<header>/) getline; }
    1 {print}
' *_${1}.txt > ${1}_combined.txt
