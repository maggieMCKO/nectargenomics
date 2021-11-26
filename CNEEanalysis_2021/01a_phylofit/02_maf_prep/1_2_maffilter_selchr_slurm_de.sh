#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-1:30:00 # 2 hours is enough for the ones with alignment
#SBATCH -p medium
#SBATCH --mem-per-cpu=5000 # doesn't need not much
#SBATCH -e logS/MF_%A.e
#SBATCH -o logS/MF_%A.o
#SBATCH -J MF

module purge
module load singularity/3.7.4 # singularity/3.2.1

# singularity pull docker://biocontainers/maffilter:v1.3.1dfsg-1b1-deb_cv1

export singularity_image="$HOME/Tools/maffilter_v1.3.1dfsg-1b1-deb_cv1.sif"
export maf_DIR="$HOME/Nectar/analysis/CNEEanalysis_2021/01a_phylofit/02_maf_prep/"
export maffilter_optionfile="${maf_DIR}running_scripts/1_1_optionfile_tmp.maffilter"

### RUN
cd ${maf_DIR}
singularity exec ${singularity_image} sh -c "cd $maf_DIR && maffilter param=${maffilter_optionfile}"
