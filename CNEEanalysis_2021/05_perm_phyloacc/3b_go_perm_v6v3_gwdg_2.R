# https://github.com/sjswuitchik/duck_comp_gen/blob/master/03_cnee_analyses/05c_go_perms.R

# v6: 
# 1. ncbi_anno > ensembl > human ensembl ortholog > go 
# 2. GRCg6a
# 3. go anno created using go_test.R  # 2. human (chicken to human orthologs to human GO) 
# v6v2
# 1. fix the ont during permutation
# the argument was 'ont' which doesn't behave as it should, ie. not really subsetting dataset, 
# so i changed the argument to 'tmp_ont'
# 2. keep GeneRatio, BgRatio for the real set
#v6v3
# use new anno
# use gene (symbol) instead of ncbi (entrez)

args = commandArgs(trailingOnly = TRUE)
nProc = args[1]

#### run GO analysis ####
# library(org.Gg.eg.db) # BiocManager::install("org.Gg.eg.db", lib)
library(tidyverse)
library(clusterProfiler)
library(parallel)

set.seed(100)

# project_path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis_2021/05_perm_phyloacc/1_rmPar/postPhyloAcc/go_perms/")
project_path = "/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/05_perm_phyloacc/1_rmPar/postPhyloAcc/go_perms/"

# custom anno
path = paste0(project_path, "/../validated.final.isoformes.galGal6.ncbi_appris_kegg_go_hsa_2022-11-21.tsv")
kegg_go_anno = read_tsv(path) # %>% 
  # mutate( gene = ifelse(gene != gene_hsa_byorg, gene_hsa_byorg, gene))
go_annoC = kegg_go_anno %>% filter(!is.na(go_id)) %>%
  dplyr::select(gene, go_id, goTerm, goTerm_def, ont) %>%
  filter(!is.na(gene)) %>%
  distinct()
dim(go_annoC) # 163449      5
length(unique(go_annoC$gene)) # 15596
sum(go_annoC$go_id =="") # 0
sum(is.na(go_annoC$go_id)) # 0

out_path = paste0(project_path, "hsa_ensembl_v6v3_2/")
dir.create(out_path)

# A. do perm ====
perm_func = function(s){
  # s = target_files[10]
  print(paste0("current: ", s))
  tmp_name = basename(s)
  tmp_name = gsub("cnee_goperms_(.*)\\.counts\\.bed", "\\1", tmp_name )
  print(paste0("current: ", tmp_name))
  
  # load in background data and clean up NCBI IDs
  bg <- read_tsv(s, col_names = F, show_col_types = FALSE) %>%
    dplyr::rename(chr = X1, start = X2, end = X3, gene = X4, combo = X5, biotype = X6, accel = X7, total = X8) %>%
    separate(combo, into = c(NA, "ncbi"), sep = "GeneID:", remove = T) %>%
    # separate(pass1, into = c("ncbi", NA), sep = "miRBase:", remove = T) %>%
    dplyr::select(-c(chr, start, end, ncbi, total, biotype))
  
  # # set up enrichment calculation
  # calc_enrich <- function(targetset, background, ont) { 
  #   enrichGO(targetset$ncbi,'org.Gg.eg.db',
  #            pvalueCutoff=1.5,
  #            qvalueCutoff = 1.5, minGSSize=5, maxGSSize=2000, 
  #            pAdjustMethod="BH",
  #            universe=background$ncbi, readable = T,
  #            keyType="ENTREZID", # keyType="SYMBOL",
  #            ont=ont) 
  # }
  
  # set up enrichment calculation
  calc_enrich <- function(targetset, background, tmp_ont) { 
    sel_custom_go_anno = go_annoC  %>% 
      filter(ont == tmp_ont)
    print(unique(sel_custom_go_anno$ont))
    term2gene = sel_custom_go_anno %>% dplyr::select(go_id, gene) %>% distinct()
    term2name = sel_custom_go_anno %>% dplyr::select(go_id, goTerm) %>% distinct()
    enricher(targetset$gene,
             pvalueCutoff=1.5,
             qvalueCutoff = 1.5, minGSSize=5, maxGSSize=2000, 
             pAdjustMethod="BH",
             universe=background$gene,
             TERM2GENE = term2gene,
             TERM2NAME = term2name) 
  }
  
  # initialize empty lists for storing permutation results
  bp.perms <- list()
  mf.perms <- list()
  cc.perms <- list()
  
  # loop through every column (each column is a permutation), select target set for given permutation, run through calc_enrich for each subontology, calculate odds ratio and enrichment scores from gene and bg ratios, and store results in list
  for (i in 1:1000) {
    print(i)
    col_name <- paste0("X", i + 8)
    perm.target <- bg %>% filter(.data[[col_name]] >= 1) %>%
      dplyr::select(gene, .data[[col_name]])
    bp <- calc_enrich(perm.target, bg, tmp_ont = "biological_process")
    mf <- calc_enrich(perm.target, bg, tmp_ont = "molecular_function")
    cc <- calc_enrich(perm.target, bg, tmp_ont = "cellular_component")
    if(! is.null(bp)){
      bp.clean <- bp@result %>% 
        separate(GeneRatio, into = c("target_in", "target_total"), sep = '/', remove= F) %>%
        separate(BgRatio, into = c("bg_in", "bg_total"), sep = '/', remove= F) %>%
        mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
               bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
               OR = target_frac/bg_frac,
               enrich = log2(target_frac/bg_frac),
               newpval = ifelse(is.na(pvalue), 1, pvalue),
               logp = -log10(newpval)) %>%
        dplyr::select(ID, GeneRatio, BgRatio, logp, target_frac, bg_frac, enrich) %>% 
        arrange(ID)
      bp.perms[[i]] <- bp.clean
    }else{bp.perms[[i]] <- NULL}
    
    if (! is.null(mf)){
      mf.clean <- mf@result %>% 
        separate(GeneRatio, into = c("target_in", "target_total"), sep = '/', remove= F) %>%
        separate(BgRatio, into = c("bg_in", "bg_total"), sep = '/', remove= F) %>%
        mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
               bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
               OR = target_frac/bg_frac,
               enrich = log2(target_frac/bg_frac),
               newpval = ifelse(is.na(pvalue), 1, pvalue),
               logp = -log10(newpval)) %>%
        dplyr::select(ID, GeneRatio, BgRatio, logp, target_frac, bg_frac, enrich) %>%
        arrange(ID)
      mf.perms[[i]] <- mf.clean
    }else{mf.perms[[i]] <- NULL}
    
    if(!is.null(cc)){
      cc.clean <- cc@result %>% 
        separate(GeneRatio, into = c("target_in", "target_total"), sep = '/', remove= F) %>%
        separate(BgRatio, into = c("bg_in", "bg_total"), sep = '/', remove= F) %>%
        mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
               bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
               OR = target_frac/bg_frac,
               enrich = log2(target_frac/bg_frac),
               newpval = ifelse(is.na(pvalue), 1, pvalue),
               logp = -log10(newpval)) %>%
        dplyr::select(ID, GeneRatio, BgRatio, logp, target_frac, bg_frac, enrich) %>%
        arrange(ID)
      cc.perms[[i]] <- cc.clean
    }else{cc.perms[[i]] <- NULL}
  }
  
  # bind permutation results
  bp.perms.clean <- bind_rows(list(bp.perms), .id = "perm")
  mf.perms.clean <- bind_rows(list(mf.perms), .id = "perm")
  cc.perms.clean <- bind_rows(list(cc.perms), .id = "perm")
  
  # write out perms so you don't have to re-run when your session inevitably crashes
  path = paste0(out_path, tmp_name, "_bp.perms.clean", Sys.Date(), ".tsv")
  write_tsv(bp.perms.clean, path, col_names = T, na = "")
  path = paste0(out_path, tmp_name, "_mf.perms.clean", Sys.Date(), ".tsv")
  write_tsv(mf.perms.clean, path, col_names = T, na = "")
  path = paste0(out_path, tmp_name, "_cc.perms.clean", Sys.Date(), ".tsv")
  write_tsv(cc.perms.clean, path, col_names = T, na = "")
  
  # filter out observed target set from background data
  target <- bg %>% filter(accel >= 1)
  
  # calculate enrichment, OR, and enrichment score for real target set for each subontology
  bp.real <- calc_enrich(target, bg, "biological_process")
  mf.real <- calc_enrich(target, bg, "molecular_function")
  cc.real <- calc_enrich(target, bg, "cellular_component")
  
  bp.real.clean <- bp.real@result %>% 
    separate(GeneRatio, into = c("target_in", "target_total"), remove = F) %>%
    separate(BgRatio, into = c("bg_in", "bg_total")) %>%
    mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
           bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
           OR = target_frac/bg_frac,
           enrich = log2(target_frac/bg_frac),
           newpval = ifelse(is.na(pvalue), 1, pvalue),
           logp = -log10(newpval)) %>%
    dplyr::select(ID, logp, target_frac, target_in, target_total, bg_frac, bg_in, bg_total, 
                  enrich, geneID) %>% # v6v2
    arrange(ID)
  mf.real.clean <- mf.real@result %>% 
    separate(GeneRatio, into = c("target_in", "target_total"), remove = F) %>%
    separate(BgRatio, into = c("bg_in", "bg_total")) %>%
    mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
           bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
           OR = target_frac/bg_frac,
           enrich = log2(target_frac/bg_frac),
           newpval = ifelse(is.na(pvalue), 1, pvalue),
           logp = -log10(newpval)) %>%
    dplyr::select(ID, logp, target_frac, target_in, target_total, bg_frac, bg_in, bg_total, 
                  enrich, geneID) %>% # v6v2
    arrange(ID)
  cc.real.clean <- cc.real@result %>% 
    separate(GeneRatio, into = c("target_in", "target_total")) %>%
    separate(BgRatio, into = c("bg_in", "bg_total")) %>%
    mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
           bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
           OR = target_frac/bg_frac,
           enrich = log2(target_frac/bg_frac),
           newpval = ifelse(is.na(pvalue), 1, pvalue),
           logp = -log10(newpval)) %>%
    dplyr::select(ID, logp, target_frac, target_in, target_total, bg_frac, bg_in, bg_total, 
                  enrich, geneID) %>% # v6v2
    arrange(ID)
  
  # merge observed results with permutation results
  if(!is_empty(bp.perms.clean) & !is_empty( bp.real.clean)){
    bp.merge <- bp.real.clean %>% 
      left_join(., bp.perms.clean, by = c("ID" = "ID"), suffix = c(".real", ".perm")) %>% 
      distinct() %>%
      arrange(ID) %>%
      dplyr::select(-perm)
    # calculate empirical pvalues (number of instances where real data is greater than or equal to permutation)
    bp.p <- bp.merge %>%
      group_by(ID) %>%
      mutate(gt_target = target_frac.perm <= target_frac.real,
             gt_enrich = enrich.perm <= enrich.real)
    # bp.cols <- sapply(bp.p[,10:11], as.numeric) # what were these 2 columns in the ori?? 0923
    bp.pval <- bp.p %>%
      group_by(ID) %>%
      mutate(sum_target = sum(gt_target) + 1,
             pVal_target = sum_target/1001,
             sum_enrich = sum(gt_enrich) + 1,
             pVal_enrich = sum_enrich/1001) %>%
      dplyr::select(ID, pVal_target, pVal_enrich, geneID, 
                    target_in, target_total, bg_in, bg_total) %>% # v6v2
      distinct()
    
    # write out GO terms
    path = paste0(out_path, tmp_name, "_bp_GOterms", Sys.Date(), ".tsv")
    write_tsv(bp.pval, path, col_names = T, na = "")
    # write out sig GO terms
    bp.pval_sig = filter(bp.pval, pVal_enrich <= 0.05) %>% 
      dplyr::select(-c(pVal_target)) %>% mutate(subontology = "biological_process") 
    path = paste0(out_path, tmp_name, "_bp_GOterms_sig", Sys.Date(), ".tsv")
    write_tsv(bp.pval_sig, path, col_names = T, na = "")
  }
  
  if(!is_empty(mf.perms.clean) & !is_empty( mf.real.clean)){
    mf.merge <- mf.real.clean %>% 
      left_join(., mf.perms.clean, by = c("ID" = "ID"), suffix = c(".real", ".perm")) %>% 
      distinct() %>%
      arrange(ID) %>%
      dplyr::select(-perm)
    # calculate empirical pvalues (number of instances where real data is greater than or equal to permutation)
    mf.p <- mf.merge %>%
      group_by(ID) %>%
      mutate(gt_target = target_frac.perm <= target_frac.real,
             gt_enrich = enrich.perm <= enrich.real)
    # mf.cols <- sapply(mf.p[,10:11], as.numeric)
    mf.pval <- mf.p %>%
      group_by(ID) %>%
      mutate(sum_target = sum(gt_target) + 1,
             pVal_target = sum_target/1001,
             sum_enrich = sum(gt_enrich) + 1,
             pVal_enrich = sum_enrich/1001) %>%
      dplyr::select(ID, pVal_target, pVal_enrich, geneID, 
                    target_in, target_total, bg_in, bg_total) %>% # v6v2
      distinct()
    
    # write out GO terms
    path = paste0(out_path, tmp_name, "_mf_GOterms", Sys.Date(), ".tsv")
    write_tsv(mf.pval, path, col_names = T, na = "")
    # write out sig GO terms
    mf.pval_sig = filter(mf.pval, pVal_enrich <= 0.05) %>% 
      dplyr::select(-c(pVal_target)) %>% 
      mutate(subontology = "molecular_function") 
    path = paste0(out_path, tmp_name, "_mf_GOterms_sig", Sys.Date(), ".tsv")
    write_tsv(mf.pval_sig, path, col_names = T, na = "")
  }
  
  if(!is_empty(cc.perms.clean) & !is_empty( cc.real.clean)){
    cc.merge <- cc.real.clean %>% 
      left_join(., cc.perms.clean, by = c("ID" = "ID"), suffix = c(".real", ".perm")) %>% 
      distinct() %>%
      arrange(ID) %>%
      dplyr::select(-perm)
    # calculate empirical pvalues (number of instances where real data is greater than or equal to permutation)
    cc.p <- cc.merge %>%
      group_by(ID) %>%
      mutate(gt_target = target_frac.perm <= target_frac.real,
             gt_enrich = enrich.perm <= enrich.real)
    # cc.cols <- sapply(cc.p[,10:11], as.numeric)
    cc.pval <- cc.p %>%
      group_by(ID) %>%
      mutate(sum_target = sum(gt_target) + 1,
             pVal_target = sum_target/1001,
             sum_enrich = sum(gt_enrich) + 1,
             pVal_enrich = sum_enrich/1001) %>%
      dplyr::select(ID, pVal_target, pVal_enrich, geneID, 
                    target_in, target_total, bg_in, bg_total) %>% # v6v2
      distinct()
    
    # write out GO terms
    path = paste0(out_path, tmp_name, "_cc_GOterms", Sys.Date(), ".tsv")
    write_tsv(cc.pval, path, col_names = T, na = "")
    # write out sig GO terms
    cc.pval_sig = filter(cc.pval, pVal_enrich <= 0.05) %>% 
      dplyr::select(-c(pVal_target)) %>% 
      mutate(subontology = "cellular_component")
    path = paste0(out_path, tmp_name, "_cc_GOterms_sig", Sys.Date(), ".tsv")
    write_tsv(cc.pval_sig, path, col_names = T, na = "")
  }
  
  
  if( exists("bp.pval_sig") ){
    sig_comb = bp.pval_sig 
    if( exists("mf.pval_sig") ){
      sig_comb = sig_comb %>% bind_rows(mf.pval_sig)
      if( exists("cc.pval_sig") ){
        sig_comb = sig_comb %>% bind_rows(cc.pval_sig)
        path = paste0(out_path, tmp_name, "_combined_GOterms_sig", Sys.Date(), ".tsv")
        write_tsv(sig_comb, path, col_names = T, na = "")
      }
    }
  }
  
  if( exists("sig_comb") ){
    return(sig_comb)
  }
  # bp.perms.clean, mf.perms.clean, cc.perms.clean
  # bp.pval, mf.pval, cc.pval,
  # bp.pval_sig, mf.pval_sig, cc.pval_sig,
  # sig_comb
}

target_files = list.files(project_path, ".counts.bed$", full.names = T); target_files
# target_files = target_files[c(4:6)]
# do_perm = lapply(target_files, perm_func)
do_perm = mclapply(X = target_files, FUN = perm_func, mc.cores = nProc) 
