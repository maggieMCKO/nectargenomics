#!/bin/bash
#SBATCH --time=00:30:00 # 1:30
#SBATCH -N 1
#SBATCH --ntasks=2
#SBATCH --mem=8G # 120
#SBATCH --partition=medium
#SBATCH -o log1_mafliftover_%A.out
#SBATCH -e log1_mafliftover_%A.err
#SBATCH -J mafliftover
##SBATCH --mail-type=ALL
##SBATCH --mail-user=mcko@orn.mpg.de
#SBATCH -C "scratch2"

source "0_setup.sh"

module purge
module load singularity/3.7.4 # singularity/3.2.1

# singularity pull docker://biocontainers/maffilter:v1.3.1dfsg-1b1-deb_cv1

export singularity_image="$HOME/Tools/maffilter_v1.3.1dfsg-1b1-deb_cv1.sif"
export maffilter_optionfile="maffilter_liftover_test/1_1_optionfile_tmp.maffilter" # direct tln
export maffilter_optionfile="maffilter_liftover_test/1_1_optionfile_tmp.maffilter2" # direct BedGraph
export maffilter_optionfile="maffilter_liftover_test/1_1_optionfile_tmp.maffilter3" # Output BedGraph

### RUN
singularity exec ${singularity_image} sh -c "cd '$HOME/Nectar/analysis/CNEEanalysis_2021/02_cnee_liftover' && maffilter param=${maffilter_optionfile}"
