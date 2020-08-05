1. prepare CNEE annotation (closest gene and cloest TSS)
    - use closest gene for GO permutation
    - run 1_cnee_anno.sh script
        - ref: https://github.com/tsackton/ratite-genomics/blob/master/04_wga/02_ce_id/postprocess_ces.sh
        - dependency:
            1. Tim's get_approx_TSS.sh
            2. setup_paths.sh: assign inputs
        - required input:
            1. ${chicken_anno} [GCF_000002315.6_GRCg6a_genomic_chrtabFixed.gff" # just chr]
                - chicken ncbi annotation, chromosome ID replaced (same as the cnee set)
            2. ${features_bed} [galGal6_final_merged_CNEEs_named_fixchr_justchr.bed" # just chr]
                - the CNEE set
        - major output:
            1. cloest_gene: galGal6_final_merged_CNEEs_cloest_genes.bed
            2. cloest_TSS: galGal6_final_merged_CNEEs_cloest_TSS.bed

2. GO-permutation
    1. prepare inputs for R script "run_enrichment_perms_gwdg.R"
        - ref:
            - https://github.com/tsackton/ratite-genomics/tree/master/07_cnee_analysis
            - https://github.com/tsackton/ratite-genomics/blob/master/07_cnee_analysis/run_enrichment_perms.R
        - dependencies:
            - library(tidyverse)
            - library(clusterProfiler) # bioconductor # BiocManager::install("clusterProfiler")
            - library(org.Gg.eg.db)    # bioconductor # BiocManager::install("org.Gg.eg.db")
            - library(parallel)
            - library(rlist)
        - required input:
            1. galGal6_final_merged_CNEEs_cloest_genes.bed
                - output of point 1
            2. phyloacc_score_postZ.tsv
                - merged output of phyloacc: merge 'score' and 'postZ' by 'No.'
                    - Rscripts/phyloacc_GO.R
        - notes:
            - you might need to edit the function 'compute_go_results' to fit the column names of your input
            - cal_enrich uses 'enrichGO' function which need ENTREZID as input
                - I use 'bitr' function to get ENTREZID from gene symbol becasue I didn't obtain ENTREZID from point 1 above.
            - 'cnee_orig_ori' is the handle for phyloacc output, you may add more conditions to restrict your interested cnees.
    2. run permutation on the cluster, using script 'run_go_perms_sbatch_m.sh'
        - ref:
            - https://github.com/tsackton/ratite-genomics/blob/master/07_cnee_analysis/run_go_perms_sbatch.sh
        - notes:
            - create an anaconda environment for running run_enrichment_perms.R with the required dependencies listed in point 2 above.
            ```bash
              module purge
              module load conda
              conda create --name R4_gwdg
              conda activate R4_gwdg
              conda install r-essentials r-base=4.0.0

              conda install -c conda-forge tidyverse

              R
              # check all requried R packages, install if lackaing
              # library(tidyverse)
              # library(clusterProfiler) # bioconductor # BiocManager::install("clusterProfiler")
              # library(org.Gg.eg.db)    # bioconductor # BiocManager::install("org.Gg.eg.db")
              # library(parallel)
              # library(rlist)
              ```

  3. compute ecdfs
    - ref:
        - https://github.com/tsackton/ratite-genomics/blob/master/07_cnee_analysis/run_gene_ecdf_calculation.R
    - run R script 'run_go_ecdf_calculation_gwdg.R'
        - in cluster
        ```bash
        Rscript run_go_ecdf_calculation_gwdg.R
        ```
        - output:
            1. goperms/original_galgal6_MF.robj
            2. goperms/original_galgal6_BP.robj

  4. analyze the output (in progress)
    - ref:
        - https://github.com/tsackton/ratite-genomics/blob/master/07_cnee_analysis/analyze_go_permutations.R
