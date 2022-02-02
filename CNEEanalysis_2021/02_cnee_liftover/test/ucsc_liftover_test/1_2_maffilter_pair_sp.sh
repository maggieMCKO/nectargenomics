#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-1:10:00 # took 50 mins
#SBATCH -p medium
#SBATCH --mem-per-cpu=8G # doesn't need not much
#SBATCH -e d1_logS/tmp_%A.e
#SBATCH -o d1_logS/tmp_%A.o
#SBATCH -J tmp

module purge
module load singularity/3.7.4 # singularity/3.2.1

# singularity pull docker://biocontainers/maffilter:v1.3.1dfsg-1b1-deb_cv1

export singularity_image="$HOME/Tools/maffilter_v1.3.1dfsg-1b1-deb_cv1.sif"
# export WD="$HOME/Nectar/analysis/CNEEanalysis_2021/02_cnee_liftover/"
export maffilter_optionfile="d1_optionfiles/1_1_optionfile_tmp.maffilter"

### RUN
# cd ${maf_DIR}
singularity exec ${singularity_image} sh -c "cd '$HOME/Nectar/analysis/CNEEanalysis_2021/02_cnee_liftover' &&
maffilter param=${maffilter_optionfile}"
