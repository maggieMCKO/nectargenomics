#!/bin/bash
#SBATCH --time=1-02:00:00 # takes 25h for nectar final run
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=4g
#SBATCH --partition=medium
#SBATCH -o log2_phylofit_%A.out
#SBATCH -e log2_phylofit_%A.err
#SBATCH -J phylofit
#SBATCH --mail-type=ALL        # BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mcko@orn.mpg.de
#SBATCH -C "scratch|scratch2"

##### SETUP INPUT and OUTPUT
export cneeAna_DIR="$HOME/Nectar/analysis/CNEEanalysis_2021/"
export input_DIR="${cneeAna_DIR}00_inputs/"
export input_tree="${input_DIR}NectarTree02.nw"

export WD="${cneeAna_DIR}01a_phylofit/03_run_phylofit/"
export tmp_DIR="${WD}01_msa/"

export concat_4ds="${WD}conca_4dSites.ss"
export neutral_model="${WD}nonconserved_4d"
export neutral_model_named="${WD}nonconserved_4d_named.mod"

##### SETUP programs
export phastPATH="$HOME/Tools/phast/bin/"
export phyloFit_PATH="${phastPATH}phyloFit"
export msa_view_PATH="${phastPATH}msa_view"
export tree_doctor_PATH="${phastPATH}tree_doctor"

##### RUN phyloFit: estimate A nonconserved phylogenetic model

${phyloFit_PATH} --tree ${input_tree} --msa-format SS --out-root ${neutral_model} ${concat_4ds}

${tree_doctor_PATH} --name-ancestors "${neutral_model}.mod" > ${neutral_model_named}
