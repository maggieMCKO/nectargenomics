#!/bin/bash
#SBATCH --time=0-00:15:00
##SBATCH -N 1 
#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=4g
#SBATCH --partition=medium
#SBATCH -o log_checkrun3_%A.o
#SBATCH -e log_check12695450_%A.e
#SBATCH -J check_csubst

# wd: /home/mpg08/mko/Nectar/analysis/csubst/01_run_csubst
export wd="$HOME/Nectar/analysis/csubst/01_run_csubst/"
export proj="out_v3"

checkFunc (){
    tmp=$1
    if ! grep -q "csubst analyze: Time elapsed" ${tmp}; then 
    	cleaned_name=$(dirname ${tmp} )
    	cleaned_name2=$(dirname ${cleaned_name} )
    	cleaned_name3=$(basename ${cleaned_name2} )
        echo -e "${cleaned_name3}"
    fi
}
export -f checkFunc
find "${wd}${proj}" -maxdepth 3 -name "log" | grep "run3" | parallel checkFunc {} > $1

# sbatch 4_check_csubst.sh unfinishedlist.txt