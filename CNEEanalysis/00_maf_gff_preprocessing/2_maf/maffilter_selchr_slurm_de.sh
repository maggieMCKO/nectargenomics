#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-2:30:00 # 2 hours is enough for the ones with alignment
#SBATCH -p medium
#SBATCH --mem-per-cpu=5000 # doesn't need not much
#SBATCH -e MF_%A.e
#SBATCH -o MF_%A.o
#SBATCH -J MF
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mcko@orn.mpg.de

module purge
module load singularity/3.2.1

# singularity pull docker://biocontainers/maffilter:v1.3.1dfsg-1b1-deb_cv1
export singularity_image="/home/mpg08/mko/Tools/maffilter_v1.3.1dfsg-1b1-deb_cv1.sif"
export maf_DIR="/home/mpg08/mko/Nectar/analysis/CNEEanalysis/00_maf_gff_preprocessing/2_maf/"
export maffilter_optionfile="${maf_DIR}optionfile_tmp.maffilter"

### RUN
cd ${maf_DIR}
singularity exec ${singularity_image} maffilter param=${maffilter_optionfile}
