#this code runs the permutations to test for GO, gene, and spatial enrichment in RARs and cRARs
library(tidyverse)
library(clusterProfiler) # bioconductor # BiocManager::install("clusterProfiler")
library(org.Gg.eg.db)    # bioconductor # BiocManager::install("org.Gg.eg.db")
# library(org.Hs.eg.db)    # bioconductor # BiocManager::install("org.Hs.eg.db")
library(parallel)
library(rlist)

#set working directory
#setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/non_coding/cnees/")

compute_go_results <- function(DF, outname, CORES, PERMS) {
  
  
  # DF = cnee_orig
  # outname = paste0(project_path, "/goperms/original_GO_run", 1)
  # CORES = 1
  # PERMS = 4
  dir.create(dirname(outname))
  
  ##INTERNAL FUNCTIONS##
  calc_enrich <- function(targetset, background, ont) { 
    enrichGO(targetset$ENTREZID,'org.Gg.eg.db',
             pvalueCutoff=1.5,
             qvalueCutoff = 1.5, minGSSize=5, maxGSSize=2000, 
             pAdjustMethod="none",
             universe=background,
             keyType="ENTREZID",
             ont=ont) 
  }
  
  get_go_perm <- function(DF, samples, golist, ont) {
    rand <- DF %>% 
      sample_n(samples) %>% 
      filter(gene != ".") %>% 
      dplyr::select(gene) %>% 
      # separate(gene, into=c("useless", "sym"), sep=":") %>% 
      separate(gene, into=c("useless", "sym"), sep="-") %>%
      distinct(sym)
    rand_ncbi = bitr(rand$sym, fromType = "SYMBOL", toType = "ENTREZID",
                           OrgDb = "org.Gg.eg.db")
    rand = rand %>% left_join(rand_ncbi, by = c("sym" = "SYMBOL"))
    
    background <- DF %>% 
      filter(gene != ".") %>% 
      dplyr::select(gene) %>% 
      # separate(gene, into=c("useless", "sym"), sep=":") %>% 
      separate(gene, into=c("useless", "sym"), sep="-") %>%
      distinct(sym)
    background_ncbi = bitr(background$sym, fromType = "SYMBOL", toType = "ENTREZID",
                     OrgDb = "org.Gg.eg.db")
    background = background %>% left_join(background_ncbi, by = c("sym" = "SYMBOL"))
    
    rand_go <- calc_enrich(targetset=rand, background=background$ENTREZID, ont=ont)
    
    
    golist %>% 
      left_join(rand_go@result, by=c("ID" = "ID")) %>% 
      separate(GeneRatio, into=c("target_in", "target_total")) %>% 
      separate(BgRatio, into=c("bg_in", "bg_total")) %>%
      mutate(newpval = ifelse(is.na(pvalue), 1, pvalue), 
             logp.perm = -log10(newpval),
             target_frac = as.numeric(target_in)/as.numeric(target_total), 
             bg_frac = as.numeric(bg_in)/as.numeric(bg_total)) %>%
      dplyr::select(ID, logp.perm, target_frac, bg_frac) %>% 
      arrange(ID)
  }
  
  #to do one permutation of the full set
  get_one_perm_set <- function(perm, input, DF, golist, ont) {
    try(lapply(input, get_go_perm, DF=DF, golist=golist, ont=ont) %>%
          dplyr::bind_rows(.id="set"), TRUE)
  }
  
  bp_perm_all <- list()
  bp_res_all <- list()
  
  mf_perm_all <- list()
  mf_res_all <- list()
  
  for (ver in c("basic")) {
    # ver = "basic"
    # c("gain", "gain_gap")
    
    cnee <- DF %>% filter(version == ver)
    
    background <- cnee %>% 
      filter(gene != ".") %>% 
      dplyr::select(gene) %>% 
      # separate(gene, into=c("ncbi", "sym"), sep=":") %>% 
      # distinct(ncbi)
      separate(gene, into=c("useless", "sym"), sep="-") %>%
      distinct(sym)
    background_ncbi = bitr(background$sym, fromType = "SYMBOL", toType = "ENTREZID",
                           OrgDb = "org.Gg.eg.db")
    background = background %>% left_join(background_ncbi, by = c("sym" = "SYMBOL"))
    
    set1 <- cnee %>% filter(gene != ".", rar) %>% 
      dplyr::select(gene) %>% 
      # separate(gene, into=c("ncbi", "sym"), sep=":") %>% 
      separate(gene, into=c("useless", "sym"), sep="-") %>%
      distinct(sym)
    set1_ncbi = bitr(set1$sym, fromType = "SYMBOL", toType = "ENTREZID",
                           OrgDb = "org.Gg.eg.db")
    set1 = set1 %>% left_join(set1_ncbi, by = c("sym" = "SYMBOL"))
    
    # set2 <- cnee %>% filter(gene != ".", crar) %>% 
    #   dplyr::select(gene) %>% 
    #   separate(gene, into=c("ncbi", "sym"), sep=":") %>% 
    #   distinct(ncbi)
    # 
    # set3 <- cnee %>% filter(gene != ".", crar_dollo) %>% 
    #   dplyr::select(gene) %>% 
    #   separate(gene, into=c("ncbi", "sym"), sep=":") %>% 
    #   distinct(ncbi)
    
    # inputs <- list("rar" = set1, "crar" = set2, "crar_dollo" = set3)
    inputs <- list("rar" = set1)
    
    bp_res_all[[ver]] <- lapply(inputs, calc_enrich, background=background$ENTREZID, ont="BP") %>% 
      lapply(slot, name="result") %>% 
      dplyr::bind_rows(.id = "set")
    
    mf_res_all[[ver]] <-  lapply(inputs, calc_enrich, background=background$ENTREZID, ont="MF") %>% 
      lapply(slot, name="result") %>% 
      dplyr::bind_rows(.id = "set")
    
    merged_mf_terms <- mf_res_all[[ver]] %>% dplyr::distinct(ID)
    merged_bp_terms <- bp_res_all[[ver]] %>% dplyr::distinct(ID)
    
    input_counts<-list("rar" = cnee %>% filter(gene != ".", rar) %>% count %>% pull(n))
                       # "crar" = cnee %>% filter(gene != ".", crar) %>% count %>% pull(n),
                       # "crar_dollo" = cnee %>% filter(gene != ".", crar_dollo) %>% count %>% pull(n))
    
    bp_perm_all[[ver]] <- mclapply(1:PERMS, get_one_perm_set, input=input_counts, DF=cnee, golist=merged_bp_terms, ont="BP", mc.cores=CORES, mc.preschedule = FALSE) %>%
      list.filter(class(.) == "data.frame") %>% 
      dplyr::bind_rows(.id="perm")
  
    mf_perm_all[[ver]] <- mclapply(1:PERMS, get_one_perm_set, input=input_counts, DF=cnee, golist=merged_mf_terms, ont="MF", mc.cores=CORES, mc.preschedule = FALSE) %>%
      list.filter(class(.) == "data.frame") %>% 
      dplyr::bind_rows(.id="perm")
    
  }
  
  bind_rows(bp_perm_all, .id="version") %>% write_tsv(paste0(outname, "_BP_perm.tsv"))
  bind_rows(bp_res_all, .id="version") %>% write_tsv(paste0(outname, "_BP_real.tsv"))              
  bind_rows(mf_perm_all, .id="version") %>% write_tsv(paste0(outname, "_MF_perm.tsv"))
  bind_rows(mf_res_all, .id="version") %>% write_tsv(paste0(outname, "_MF_real.tsv"))       
}

### REAL WORK ###

args <- commandArgs(trailingOnly = TRUE)
#args are 1 annotation file name, 2 permutation index ID for slurm batch processing, 3 number of cores, 4 number of permutations, 5 is data path

path_to_data <- args[5]

# project_path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/2_GO/0_prepare/")
project_path = paste0("/home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/2_GO/0_prepare/") # gwdg

path = paste0(project_path, "galGal6_final_merged_CNEEs_cloest_genes.bed")
# gene_gg = read_tsv(paste0("../04_wga/03_ce_annotation/cnees.", args[1], ".annotation"), col_names = c("cnee", "gene"))
gene_gg = read_tsv(path, col_names = F)
gene_gg = gene_gg[, c(4, 8)]
names(gene_gg) = c("cnee", "gene")
# there are duplicates in gene_gg

# path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/1_runphyloacc/1_run_phyloacc/phyloacc_score_postZ.tsv")
path = paste0("/home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/1_run_phyloacc/phyloacc_score_postZ.tsv") # gwdg
cnee_orig_ori <- read_tsv(path) 
names(cnee_orig_ori) = gsub("ID", "cnee", names(cnee_orig_ori))

cnee_orig = cnee_orig_ori %>% 
  dplyr::select(No., cnee, logBF1, logBF2) %>%
  full_join(gene_gg, by=c("cnee" = "cnee")) %>%
  mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 , TRUE, FALSE)) %>%
  distinct(cnee, .keep_all=TRUE) %>% 
  dplyr::select(cnee, rar, gene) %>%
  mutate(version = "basic")
rm(cnee_orig_ori)

compute_go_results(cnee_orig, paste0(path_to_data, "/goperms/original_GO_", args[1], "_run", args[2]), args[3], args[4])

# test
# project_path = paste0("/home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/2_GO/1_enrichment/") # gwdg
# compute_go_results(cnee_orig, paste0(project_path, "/goperms/original_GO_run", 1), 1, 4)

# cnee_orig <- read_tsv(paste0(path_to_data, "/final_original_cnee.tsv.gz")) %>% 
#   dplyr::select(version, cnee, logBF1, logBF2, it_pp_loss, ti_pp_loss, neo_tip_loss, floss_cl_pp, floss_cl_pp_dollo) %>%
#   full_join(gene_gg, by=c("cnee" = "cnee")) %>%
#   mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 & (it_pp_loss + ti_pp_loss) < 1 & neo_tip_loss < 1, TRUE, FALSE),
#          crar = ifelse(rar & floss_cl_pp >= 1.8, TRUE, FALSE),
#          crar_dollo = ifelse(rar & floss_cl_pp_dollo >= 1.8, TRUE, FALSE)) %>%
#   distinct(cnee, version, .keep_all=TRUE) %>%
#   dplyr::select(cnee, version, rar, crar, crar_dollo, gene)

# compute_go_results(cnee_orig, paste0(path_to_data, "/goperms/original_GO_", args[1], "_run", args[2]), args[3], args[4])
