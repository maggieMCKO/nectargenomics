#!/bin/bash
#SBATCH --time=9:00:00
#SBATCH -c 1
##SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4g
#SBATCH --partition=medium
#SBATCH -o phylop_%A.out
#SBATCH -e phylop_%A.err
#SBATCH -J phylop
#SBATCH --mail-type=ALL        # BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mcko@orn.mpg.de

##### SETUP INPUT and OUTPUT
cneeAna_DIR="/home/mpg08/mko/Nectar/analysis/CNEEanalysis/"

maf_DIR="${cneeAna_DIR}MAF/"
input_maf="${maf_DIR}multi.anno.maf"

phylofit_DIR="${cneeAna_DIR}01_phylofit/01_msa3/"
neutral_model="${phylofit_DIR}nonconserved_4d.mod"

phylop_DIR="${cneeAna_DIR}02_phylop/"
features_bed="${phylop_DIR}galGal6_final_merged_CNEEs_named_fixchr2.bed"
features_out="${phylop_DIR}features.out"

WD="${cneeAna_DIR}02_phylop/subtree/"
neutral_model_named="${WD}nonconserved_4d_named.mod"

##### SETUP programs
export phastPATH="/home/mpg08/mko/Tools/phast/bin/"
export phyloP_PATH="${phastPATH}phyloP"
export tree_doctor_PATH="${phastPATH}tree_doctor"


##### RUN
### Run phyloP with the subtree and branch options.
# ${tree_doctor_PATH} --name-ancestors ${neutral_model} > ${neutral_model_named}

# Scores describing lineage specific conservation can then be computed using the --subtree option. Here we use the --base-by-base option which outputs multiple values per site, in a method-dependent way.
${phyloP_PATH} --method LRT --subtree HLcalAnn5-HLfloFus1 --mode CONACC --msa-format MAF --base-by-base ${neutral_model_named} ${input_maf} > subtree_Hummingbird_basebybase

# Similarly a features file can be provide along with the subtree option to get a table of p-values and related statistics with one row per feature.
${phyloP_PATH} --method LRT --subtree HLcalAnn5-HLfloFus1 --mode CONACC --features ${features_bed} --msa-format MAF ${neutral_model_named} ${input_maf} > subtree_Hummingbird_basebybase_cnee


# The --branch option is similar to --subtree, but it partitions the tree into the set of named branches, and all the remaining branches before testing for conservation/acceleration in the set of names branches relative to the others. A comma-delimited list of child nodes can be provided as an argument.

${phyloP_PATH}  --method LRT --branch HLcalAnn5-HLfloFus1 --mode CONACC --wig-scores --msa-format MAF ${neutral_model_named} ${input_maf} > branch_Hummingbird.wig

${phyloP_PATH}  --method LRT --branch HLcalAnn5-HLfloFus1 --mode CONACC --wig-scores --msa-format MAF --features ${features_bed} ${neutral_model_named} ${input_maf} > branch_Hummingbird_cnee.wig


# USAGE: phyloP [OPTIONS] tree.mod [alignment] > out
