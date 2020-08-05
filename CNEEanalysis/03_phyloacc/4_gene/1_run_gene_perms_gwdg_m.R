library(tidyverse)
library(parallel)

### GENE ENRICHMENT HERE ###
#for gene enrichment tests, the idea is to randomly permute each set and get count of CNEEs per gene
#strategy here is to make an indicator variable, compute "real" T/F per gene, and then shuffle indicator variable

get_gene_counts <- function(DF, indicator) {
  #  indcol <- enquo(indicator)
  #  indcol <- indicator
  DF %>% filter(gene != ".") %>% mutate(in_target = !!indicator) %>% 
    count(gene, in_target) %>% filter(!is.na(in_target)) %>% spread(in_target, n, fill=0, drop=FALSE, sep="_")
}

perm_gene_counts <- function(perm, DF, indicator) {
  #  indcol <- enquo(indicator)
  DF %>% filter(gene != ".") %>% mutate(rand = sample(!!indicator)) %>% 
    count(gene, rand) %>% filter(!is.na(rand)) %>% spread(rand, n, fill=0, drop=FALSE, sep="_")
}

#function to do work

compute_gene_results <- function(DF, outname, CORES, PERMS) {
  # DF = cnee_orig
  # outname = paste0(project_path, "/goperms/original_GO_run", 1)
  # CORES = 1
  # PERMS = 4
  dir.create(dirname(outname))
  
  gene_perm_all <- list()
  gene_res_all <- list()
  
  for (ver in c("basic")) {
    
    cnee <- DF %>% filter(version == ver)
    gene_res_all[[ver]] <- lapply(c(quo(rar), quo(filter_non_target_5sp), quo(filter_non_target_humm2anc), 
                                    quo(filter_non_target_humm2anc_honeyanc)), 
                                  get_gene_counts, DF=cnee) %>% bind_rows(.id="set")
  
    para_cores <- CORES
    num_perms <- PERMS
  
    gene_counts_perm_list <- list(rar = mclapply(1:num_perms, perm_gene_counts, DF=cnee, indicator=quo(rar), 
                                                 mc.preschedule = FALSE, mc.cores = para_cores) %>% bind_rows(.id="perm"),
                                  # 1
                                  filter_non_target_5sp = mclapply(1:num_perms, perm_gene_counts, DF=cnee, 
                                                                   indicator=quo(filter_non_target_5sp), 
                                                                   mc.preschedule = FALSE, mc.cores = para_cores) %>% bind_rows(.id="perm"),
                                  # 2
                                  filter_non_target_humm2anc = mclapply(1:num_perms, perm_gene_counts, DF=cnee, 
                                                                        indicator=quo(filter_non_target_humm2anc), 
                                                                        mc.preschedule = FALSE, mc.cores = para_cores) %>% bind_rows(.id="perm"),
                                  # # 2.2
                                  # filter_non_target_humm2 = mclapply(1:num_perms, perm_gene_counts, DF=cnee, 
                                  #                                  indicator=quo(filter_non_target_humm2), 
                                  #                                  mc.preschedule = FALSE, mc.cores = para_cores) %>% bind_rows(.id="perm"),
                                  # # 3
                                  # filter_non_target_humm2anc_triMol = mclapply(1:num_perms, perm_gene_counts, DF=cnee, 
                                  #                                       indicator=quo(filter_non_target_humm2anc_triMol), 
                                  #                                       mc.preschedule = FALSE, mc.cores = para_cores) %>% bind_rows(.id="perm"),
                                  # 4
                                  filter_non_target_humm2anc_honeyanc = mclapply(1:num_perms, perm_gene_counts, DF=cnee,
                                                                   indicator=quo(filter_non_target_humm2anc_honeyanc),
                                                                   mc.preschedule = FALSE, mc.cores = para_cores) %>% bind_rows(.id="perm"))
                                  # # 5
                                  # filter_non_target_honey2anc_triMol = mclapply(1:num_perms, perm_gene_counts, DF=cnee, 
                                  #                                       indicator=quo(filter_non_target_honey2anc_triMol), 
                                  #                                       mc.preschedule = FALSE, mc.cores = para_cores) %>% bind_rows(.id="perm"))
    
    gene_perm_all[[ver]] <- gene_counts_perm_list %>% bind_rows(.id="set")
  
  }

  bind_rows(gene_perm_all, .id="version") %>% write_tsv(paste0(outname, "_perm.tsv"))
  bind_rows(gene_res_all, .id="version") %>% write_tsv(paste0(outname, "_real.tsv"))
                                                                                                              
}

args <- commandArgs(trailingOnly = TRUE)
#args are 1 annotation file name, 2 permutation index ID for slurm batch processing, 3 number of cores, 4 number of permutations, 5 is data path

path_to_data <- args[5]

# project_path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/2_GO/0_prepare/")
project_path = paste0("/home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/2_GO/0_prepare/") # gwdg

path = paste0(project_path, "galGal6_final_merged_CNEEs_cloest_genes.bed")
gene_gg = read_tsv(path, col_names = F)
gene_gg = gene_gg[, c(4, 8)]
names(gene_gg) = c("cnee", "gene")
# there are duplicates in gene_gg

# path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/1_runphyloacc/1_run_phyloacc/phyloacc_score_postZ.tsv")
path = paste0("/home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/1_run_phyloacc/phyloacc_score_postZ.tsv.gz") # gwdg
cnee_orig_ori <- read_tsv(path) 
names(cnee_orig_ori) = gsub("ID", "cnee", names(cnee_orig_ori))
names(cnee_orig_ori)[250:257]

cnee_orig = cnee_orig_ori %>% 
  dplyr::select(No., cnee, logBF1, logBF2,
                filter_non_target_5sp, filter_non_target_humm2anc, filter_non_target_humm2, 
                filter_non_target_humm2anc_triMol, filter_non_target_humm2anc_honeyanc, filter_non_target_honey2anc_triMol) %>%
  full_join(gene_gg, by=c("cnee" = "cnee")) %>%
  mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 , TRUE, FALSE)) %>%
  distinct(cnee, .keep_all=TRUE) %>% 
  dplyr::select(cnee, gene, rar,
  filter_non_target_5sp, filter_non_target_humm2anc, filter_non_target_humm2, 
  filter_non_target_humm2anc_triMol, filter_non_target_humm2anc_honeyanc, filter_non_target_honey2anc_triMol) %>%
  mutate(version = "basic")
rm(cnee_orig_ori)

# cnee_orig <- read_tsv(paste0(path_to_data, "/final_original_cnee.tsv.gz")) %>% 
#   select(version, cnee, logBF1, logBF2, it_pp_loss, ti_pp_loss, neo_tip_loss, floss_cl_pp, floss_cl_pp_dollo) %>%
#   full_join(gene_gg, by=c("cnee" = "cnee")) %>%
#   mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 & (it_pp_loss + ti_pp_loss) < 1 & neo_tip_loss < 1, TRUE, FALSE),
#          crar = ifelse(rar & floss_cl_pp >= 1.8, TRUE, FALSE),
#          crar_dollo = ifelse(rar & floss_cl_pp_dollo >= 1.8, TRUE, FALSE)) %>%
#   distinct(cnee, version, .keep_all=TRUE) %>%
#   select(cnee, version, rar, crar, crar_dollo, gene)


#note in ext2, convergence defined as ratites + cormorants
#cnee_ext2 <- read_tsv("final_extended_cnee.tsv.gz") %>% 
#  select(version, cnee, logBF1, logBF2, it_pp_loss, ti_pp_loss, neo_tip_loss, floss_cl_pp, floss_cl_pp_dollo, gc_pp_loss) %>%
#  full_join(gene_gg, by=c("cnee" = "cnee")) %>%
#  mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 & (it_pp_loss + ti_pp_loss) < 1 & neo_tip_loss < 1, TRUE, FALSE),
#         crar = ifelse(rar & floss_cl_pp >= 1.8 & gc_pp_loss > 0.90, TRUE, FALSE),
#         crar_dollo = ifelse(rar & floss_cl_pp_dollo >= 1.8 & gc_pp_loss > 0.90, TRUE, FALSE)) %>%
#  distinct(cnee, version, .keep_all=TRUE) %>%
#  select(cnee, version, rar, crar, crar_dollo, gene)

compute_gene_results(cnee_orig, paste0(path_to_data, "/geneperms/original_gene_", args[1], "_run", args[2]), args[3], args[4])
# compute_gene_results(cnee_ext, paste0(path_to_data, "/geneperms/extended_gene_", args[1], "_run", args[2]), args[3], args[4])
# #compute_gene_results(cnee_ext2, paste0("geneperms/extended_ratiteVcorm_gene_", args[1], "_run", args[2]), args[3], args[4])
# compute_gene_results(cnee_red, paste0(path_to_data, "/geneperms/reduced_gene_", args[1], "_run", args[2]), args[3], args[4])


