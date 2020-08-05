1. convergence-permutation
    1. prepare inputs for R script "run_enrichment_perms_gwdg.R"
        - ref:
            - https://github.com/tsackton/ratite-genomics/tree/master/07_cnee_analysis
            - https://github.com/tsackton/ratite-genomics/blob/master/07_cnee_analysis/run_convergence_perms.R
        - dependencies:
            - library(tidyverse)
            - library(parallel)
            - library(ape)
        - required input:
            1. phyloacc_score_postZ.tsv
                - merged output of phyloacc: merge 'score' and 'postZ' by 'No.'
                - only 'acc' state is needed
        - notes:
            -
    2. run permutation on the cluster, using script 'run_convergence_perms_sbatch_m.sh'
        - ref:
            - https://github.com/tsackton/ratite-genomics/blob/master/07_cnee_analysis/run_convergence_perms_sbatch.sh
        - notes:
            - create an anaconda environment for running 'run_convergence_perms_sbatch_m.sh' with the required dependencies listed in point 2 above.
            ```bash
                # need update, see condaR.md
              module purge
              module load conda
              conda create --name R4_gwdg
              conda activate R4_gwdg
              conda install r-essentials r-base=4.0.0

              conda install -c conda-forge tidyverse

              R
              # check all requried R packages, install if lackaing
              # library(tidyverse)
              # library(parallel)
              # library(rlist)
              ```

  3. analyze the output 'analyze_convergence_permutations_m.R'
    - ref:
        - https://github.com/tsackton/ratite-genomics/blob/master/07_cnee_analysis/analyze_convergence_permutations.R
    - part 1, i did it in the cluster
