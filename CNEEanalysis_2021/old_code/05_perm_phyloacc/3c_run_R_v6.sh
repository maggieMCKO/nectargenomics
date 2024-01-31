#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=4g
#SBATCH --partition=medium
#SBATCH -o log3c_perm_go_R_v6_%A.o
#SBATCH -e log3c_perm_go_R_v6_%A.e
#SBATCH -J perm_go_Rv6

# WD: /home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/05_perm_phyloacc/1_rmPar


module purge
module load anaconda2/2019.10
source activate R_phylo

# Rscript 3b_go_perm_v6v2_gwdg.R 16
# Rscript 3b_go_perm_v6v3_gwdg.R 20
Rscript 3b_go_perm_v6v3_gwdg_2.R 20
# Rscript 3b_go_perm_v6v3_gwdg_2_cladespecific.R 20
