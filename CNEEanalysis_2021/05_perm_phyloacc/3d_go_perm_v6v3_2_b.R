# https://github.com/sjswuitchik/duck_comp_gen/blob/master/03_cnee_analyses/05c_go_perms.R

# v6: 
# 1. ncbi_anno > ensembl > human ensembl ortholog > go 
# 2. GRCg6a
# 3. go anno created using go_test.R  # 2. human (chicken to human orthologs to human GO) 
# v6v2
# 1. fix the ont during permutation
# the argument was 'ont' which doesn't behave as it should, ie. not really subsetting dataset, 
# so i changed the argument to 'tmp_ont'
# 2. keep GeneRatio (target_in, target_total), BgRatio (bg_in, bg_total) for the real set
# (they were kept before, just i didn't realize it)
# 3. it seems like the p values we used until now were not adjusted for multiple comparison
#v6v3
# use new anno
# use gene (symbol) instead of ncbi (entrez)
# v6v3_2
# use new anno 1118 version

#### run GO analysis ####
library(tidyverse)
library(clusterProfiler)
library(RColorBrewer)
library(parallel)
library(enrichplot)

openfile <- function(filepath){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  type = tolower(os);
  switch(os, Linux = system(paste0("xdg-open ", filepath)),
         Windows = system(paste0("open \"", filepath, "\"")),
         osx = system(paste0("open \"", filepath, "\""))
  )}

FontSize = 2.5
# AxisTxFontSizeSize = 8
# AxisTxFontSizeSize_s = 6
# AxisTitleFontSizeSize = 10
AxisTxFontSizeSize = 6
AxisTxFontSizeSize_s = 5
AxisTitleFontSizeSize = 6
Factor_mmtoin = 0.0393701
Width_HalfCol = 85*Factor_mmtoin 

theme_m = theme(panel.background = element_blank(),
                plot.background = element_blank(), # defualt is white
                plot.title = element_text(colour = "black", size = AxisTitleFontSizeSize, hjust = 0.5),
                plot.margin = margin(0.5, 2, 0.5, 0.5, "line"),
                panel.grid = element_blank(),
                # panel.grid.major.x = element_line(color = "grey", size = .2),
                # panel.grid.minor.x = element_blank(),
                # panel.grid.major.y = element_line(color = "grey", size = .2),
                # panel.grid.minor.y = element_blank(), #element_line(color = "red", size = .2), 
                panel.border = element_rect(color = "black", fill= NA, size =.2), 
                # axis.line = element_line(color = "black", size = 0.5),
                axis.title.x = element_text(colour = "black", size = AxisTitleFontSizeSize, hjust = 0.5),
                axis.title.y = element_text(colour = "black", size = AxisTitleFontSizeSize),
                axis.text.x = element_text(colour = "black", size = AxisTxFontSizeSize),
                axis.text.y = element_text(colour = "black", size = AxisTxFontSizeSize),
                axis.ticks = element_line(colour = "black", size = 0.05),
                strip.text = element_text(colour = "black", size = AxisTxFontSizeSize),
                legend.background = element_blank(),
                legend.key = element_blank(),
                legend.key.size = unit(8, "pt"),
                legend.title = element_blank(),
                legend.text = element_text(colour = "black", size = AxisTxFontSizeSize),
                legend.position = "right")
set.seed(100)

project_path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis_2021/05_perm_phyloacc/1_rmPar/postPhyloAcc/go_perms/")

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

# # A. do perm ====
# perm_func = function(s){
#   # s = target_files[10]
#   print(paste0("current: ", s))
#   tmp_name = basename(s)
#   tmp_name = gsub("cnee_goperms_(.*)\\.counts\\.bed", "\\1", tmp_name )
#   print(paste0("current: ", tmp_name))
#   
#   # load in background data and clean up NCBI IDs
#   bg <- read_tsv(s, col_names = F, show_col_types = FALSE) %>%
#     dplyr::rename(chr = X1, start = X2, end = X3, gene = X4, combo = X5, biotype = X6, accel = X7, total = X8) %>%
#     separate(combo, into = c(NA, "ncbi"), sep = "GeneID:", remove = T) %>%
#     # separate(pass1, into = c("ncbi", NA), sep = "miRBase:", remove = T) %>%
#     dplyr::select(-c(chr, start, end, ncbi, total, biotype))
#   
#   # # set up enrichment calculation
#   # calc_enrich <- function(targetset, background, ont) { 
#   #   enrichGO(targetset$ncbi,'org.Gg.eg.db',
#   #            pvalueCutoff=1.5,
#   #            qvalueCutoff = 1.5, minGSSize=5, maxGSSize=2000, 
#   #            pAdjustMethod="BH",
#   #            universe=background$ncbi, readable = T,
#   #            keyType="ENTREZID", # keyType="SYMBOL",
#   #            ont=ont) 
#   # }
#   
#   # set up enrichment calculation
#   calc_enrich <- function(targetset, background, tmp_ont) { 
#     sel_custom_go_anno = go_annoC  %>% 
#       filter(ont == tmp_ont)
#     print(unique(sel_custom_go_anno$ont))
#     term2gene = sel_custom_go_anno %>% dplyr::select(go_id, gene) %>% distinct()
#     term2name = sel_custom_go_anno %>% dplyr::select(go_id, goTerm) %>% distinct()
#     enricher(targetset$gene,
#              pvalueCutoff=1.5,
#              qvalueCutoff = 1.5, minGSSize=5, maxGSSize=2000, 
#              pAdjustMethod="BH",
#              universe=background$gene,
#              TERM2GENE = term2gene,
#              TERM2NAME = term2name) 
#   }
#   
#   # initialize empty lists for storing permutation results
#   bp.perms <- list()
#   mf.perms <- list()
#   cc.perms <- list()
#   
#   # loop through every column (each column is a permutation), select target set for given permutation, run through calc_enrich for each subontology, calculate odds ratio and enrichment scores from gene and bg ratios, and store results in list
#   for (i in 1:1000) {
#     print(i)
#     col_name <- paste0("X", i + 8)
#     perm.target <- bg %>% filter(.data[[col_name]] >= 1) %>%
#       dplyr::select(gene, .data[[col_name]])
#     bp <- calc_enrich(perm.target, bg, tmp_ont = "biological_process")
#     mf <- calc_enrich(perm.target, bg, tmp_ont = "molecular_function")
#     cc <- calc_enrich(perm.target, bg, tmp_ont = "cellular_component")
#     if(! is.null(bp)){
#       bp.clean <- bp@result %>% 
#         separate(GeneRatio, into = c("target_in", "target_total"), sep = '/', remove= F) %>%
#         separate(BgRatio, into = c("bg_in", "bg_total"), sep = '/', remove= F) %>%
#         mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
#                bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
#                OR = target_frac/bg_frac,
#                enrich = log2(target_frac/bg_frac),
#                newpval = ifelse(is.na(pvalue), 1, pvalue),
#                logp = -log10(newpval)) %>%
#         dplyr::select(ID, GeneRatio, BgRatio, logp, target_frac, bg_frac, enrich) %>% 
#         arrange(ID)
#       bp.perms[[i]] <- bp.clean
#     }else{bp.perms[[i]] <- NULL}
#     
#     if (! is.null(mf)){
#       mf.clean <- mf@result %>% 
#         separate(GeneRatio, into = c("target_in", "target_total"), sep = '/', remove= F) %>%
#         separate(BgRatio, into = c("bg_in", "bg_total"), sep = '/', remove= F) %>%
#         mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
#                bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
#                OR = target_frac/bg_frac,
#                enrich = log2(target_frac/bg_frac),
#                newpval = ifelse(is.na(pvalue), 1, pvalue),
#                logp = -log10(newpval)) %>%
#         dplyr::select(ID, GeneRatio, BgRatio, logp, target_frac, bg_frac, enrich) %>%
#         arrange(ID)
#       mf.perms[[i]] <- mf.clean
#     }else{mf.perms[[i]] <- NULL}
#     
#     if(!is.null(cc)){
#       cc.clean <- cc@result %>% 
#         separate(GeneRatio, into = c("target_in", "target_total"), sep = '/', remove= F) %>%
#         separate(BgRatio, into = c("bg_in", "bg_total"), sep = '/', remove= F) %>%
#         mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
#                bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
#                OR = target_frac/bg_frac,
#                enrich = log2(target_frac/bg_frac),
#                newpval = ifelse(is.na(pvalue), 1, pvalue),
#                logp = -log10(newpval)) %>%
#         dplyr::select(ID, GeneRatio, BgRatio, logp, target_frac, bg_frac, enrich) %>%
#         arrange(ID)
#       cc.perms[[i]] <- cc.clean
#     }else{cc.perms[[i]] <- NULL}
#   }
#   
#   # bind permutation results
#   bp.perms.clean <- bind_rows(list(bp.perms), .id = "perm")
#   mf.perms.clean <- bind_rows(list(mf.perms), .id = "perm")
#   cc.perms.clean <- bind_rows(list(cc.perms), .id = "perm")
#   
#   # write out perms so you don't have to re-run when your session inevitably crashes
#   path = paste0(out_path, tmp_name, "_bp.perms.clean", Sys.Date(), ".tsv")
#   write_tsv(bp.perms.clean, path, col_names = T, na = "")
#   path = paste0(out_path, tmp_name, "_mf.perms.clean", Sys.Date(), ".tsv")
#   write_tsv(mf.perms.clean, path, col_names = T, na = "")
#   path = paste0(out_path, tmp_name, "_cc.perms.clean", Sys.Date(), ".tsv")
#   write_tsv(cc.perms.clean, path, col_names = T, na = "")
#   
#   # filter out observed target set from background data
#   target <- bg %>% filter(accel >= 1)
#   
#   # calculate enrichment, OR, and enrichment score for real target set for each subontology
#   bp.real <- calc_enrich(target, bg, "biological_process")
#   mf.real <- calc_enrich(target, bg, "molecular_function")
#   cc.real <- calc_enrich(target, bg, "cellular_component")
#   
#   bp.real.clean <- bp.real@result %>% 
#     separate(GeneRatio, into = c("target_in", "target_total"), remove = F) %>%
#     separate(BgRatio, into = c("bg_in", "bg_total")) %>%
#     mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
#            bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
#            OR = target_frac/bg_frac,
#            enrich = log2(target_frac/bg_frac),
#            newpval = ifelse(is.na(pvalue), 1, pvalue),
#            logp = -log10(newpval)) %>%
#     dplyr::select(ID, logp, target_frac, target_in, target_total, bg_frac, bg_in, bg_total, 
#                   enrich, geneID) %>% # v6v2
#     arrange(ID)
#   mf.real.clean <- mf.real@result %>% 
#     separate(GeneRatio, into = c("target_in", "target_total"), remove = F) %>%
#     separate(BgRatio, into = c("bg_in", "bg_total")) %>%
#     mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
#            bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
#            OR = target_frac/bg_frac,
#            enrich = log2(target_frac/bg_frac),
#            newpval = ifelse(is.na(pvalue), 1, pvalue),
#            logp = -log10(newpval)) %>%
#     dplyr::select(ID, logp, target_frac, target_in, target_total, bg_frac, bg_in, bg_total, 
#                   enrich, geneID) %>% # v6v2
#     arrange(ID)
#   cc.real.clean <- cc.real@result %>% 
#     separate(GeneRatio, into = c("target_in", "target_total")) %>%
#     separate(BgRatio, into = c("bg_in", "bg_total")) %>%
#     mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
#            bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
#            OR = target_frac/bg_frac,
#            enrich = log2(target_frac/bg_frac),
#            newpval = ifelse(is.na(pvalue), 1, pvalue),
#            logp = -log10(newpval)) %>%
#     dplyr::select(ID, logp, target_frac, target_in, target_total, bg_frac, bg_in, bg_total, 
#                   enrich, geneID) %>% # v6v2
#     arrange(ID)
#   
#   # merge observed results with permutation results
#   if(!is_empty(bp.perms.clean) & !is_empty( bp.real.clean)){
#     bp.merge <- bp.real.clean %>% 
#       left_join(., bp.perms.clean, by = c("ID" = "ID"), suffix = c(".real", ".perm")) %>% 
#       distinct() %>%
#       arrange(ID) %>%
#       dplyr::select(-perm)
#     # calculate empirical pvalues (number of instances where real data is greater than or equal to permutation)
#     bp.p <- bp.merge %>%
#       group_by(ID) %>%
#       mutate(gt_target = target_frac.perm <= target_frac.real,
#              gt_enrich = enrich.perm <= enrich.real)
#     # bp.cols <- sapply(bp.p[,10:11], as.numeric) # what were these 2 columns in the ori?? 0923
#     bp.pval <- bp.p %>%
#       group_by(ID) %>%
#       mutate(sum_target = sum(gt_target) + 1,
#              pVal_target = sum_target/1001,
#              sum_enrich = sum(gt_enrich) + 1,
#              pVal_enrich = sum_enrich/1001) %>%
#       dplyr::select(ID, pVal_target, pVal_enrich, geneID, 
#                     target_in, target_total, bg_in, bg_total) %>% # v6v2
#       distinct()
#     
#     # write out GO terms
#     path = paste0(out_path, tmp_name, "_bp_GOterms", Sys.Date(), ".tsv")
#     write_tsv(bp.pval, path, col_names = T, na = "")
#     # write out sig GO terms
#     bp.pval_sig = filter(bp.pval, pVal_enrich <= 0.05) %>% 
#       dplyr::select(-c(pVal_target)) %>% mutate(subontology = "biological_process") 
#     path = paste0(out_path, tmp_name, "_bp_GOterms_sig", Sys.Date(), ".tsv")
#     write_tsv(bp.pval_sig, path, col_names = T, na = "")
#   }
#   
#   if(!is_empty(mf.perms.clean) & !is_empty( mf.real.clean)){
#     mf.merge <- mf.real.clean %>% 
#       left_join(., mf.perms.clean, by = c("ID" = "ID"), suffix = c(".real", ".perm")) %>% 
#       distinct() %>%
#       arrange(ID) %>%
#       dplyr::select(-perm)
#     # calculate empirical pvalues (number of instances where real data is greater than or equal to permutation)
#     mf.p <- mf.merge %>%
#       group_by(ID) %>%
#       mutate(gt_target = target_frac.perm <= target_frac.real,
#              gt_enrich = enrich.perm <= enrich.real)
#     # mf.cols <- sapply(mf.p[,10:11], as.numeric)
#     mf.pval <- mf.p %>%
#       group_by(ID) %>%
#       mutate(sum_target = sum(gt_target) + 1,
#              pVal_target = sum_target/1001,
#              sum_enrich = sum(gt_enrich) + 1,
#              pVal_enrich = sum_enrich/1001) %>%
#       dplyr::select(ID, pVal_target, pVal_enrich, geneID, 
#                     target_in, target_total, bg_in, bg_total) %>% # v6v2
#       distinct()
#     
#     # write out GO terms
#     path = paste0(out_path, tmp_name, "_mf_GOterms", Sys.Date(), ".tsv")
#     write_tsv(mf.pval, path, col_names = T, na = "")
#     # write out sig GO terms
#     mf.pval_sig = filter(mf.pval, pVal_enrich <= 0.05) %>% 
#       dplyr::select(-c(pVal_target)) %>% 
#       mutate(subontology = "molecular_function") 
#     path = paste0(out_path, tmp_name, "_mf_GOterms_sig", Sys.Date(), ".tsv")
#     write_tsv(mf.pval_sig, path, col_names = T, na = "")
#   }
#   
#   if(!is_empty(cc.perms.clean) & !is_empty( cc.real.clean)){
#     cc.merge <- cc.real.clean %>% 
#       left_join(., cc.perms.clean, by = c("ID" = "ID"), suffix = c(".real", ".perm")) %>% 
#       distinct() %>%
#       arrange(ID) %>%
#       dplyr::select(-perm)
#     # calculate empirical pvalues (number of instances where real data is greater than or equal to permutation)
#     cc.p <- cc.merge %>%
#       group_by(ID) %>%
#       mutate(gt_target = target_frac.perm <= target_frac.real,
#              gt_enrich = enrich.perm <= enrich.real)
#     # cc.cols <- sapply(cc.p[,10:11], as.numeric)
#     cc.pval <- cc.p %>%
#       group_by(ID) %>%
#       mutate(sum_target = sum(gt_target) + 1,
#              pVal_target = sum_target/1001,
#              sum_enrich = sum(gt_enrich) + 1,
#              pVal_enrich = sum_enrich/1001) %>%
#       dplyr::select(ID, pVal_target, pVal_enrich, geneID, 
#                     target_in, target_total, bg_in, bg_total) %>% # v6v2
#       distinct()
#     
#     # write out GO terms
#     path = paste0(out_path, tmp_name, "_cc_GOterms", Sys.Date(), ".tsv")
#     write_tsv(cc.pval, path, col_names = T, na = "")
#     # write out sig GO terms
#     cc.pval_sig = filter(cc.pval, pVal_enrich <= 0.05) %>% 
#       dplyr::select(-c(pVal_target)) %>% 
#       mutate(subontology = "cellular_component")
#     path = paste0(out_path, tmp_name, "_cc_GOterms_sig", Sys.Date(), ".tsv")
#     write_tsv(cc.pval_sig, path, col_names = T, na = "")
#   }
#   
#   
#   if( exists("bp.pval_sig") ){
#     sig_comb = bp.pval_sig 
#     if( exists("mf.pval_sig") ){
#       sig_comb = sig_comb %>% bind_rows(mf.pval_sig)
#       if( exists("cc.pval_sig") ){
#         sig_comb = sig_comb %>% bind_rows(cc.pval_sig)
#         path = paste0(out_path, tmp_name, "_combined_GOterms_sig", Sys.Date(), ".tsv")
#         write_tsv(sig_comb, path, col_names = T, na = "")
#       }
#     }
#   }
#   
#   if( exists("sig_comb") ){
#     return(sig_comb)
#   }
#   # bp.perms.clean, mf.perms.clean, cc.perms.clean
#   # bp.pval, mf.pval, cc.pval,
#   # bp.pval_sig, mf.pval_sig, cc.pval_sig,
#   # sig_comb
# }
# 
# target_files = list.files(project_path, ".counts.bed$", full.names = T); target_files
# # target_files = target_files[c(4:6)]
# # do_perm = lapply(target_files, perm_func)
# do_perm = mclapply(X = target_files, FUN = perm_func, mc.cores = nProc) 

# B. analysis - sig ####
library(biomaRt) # installed dev version with BiocManager::install('grimbough/biomaRt')
library(grid)
library(gridExtra)

sig_files = list.files(out_path, "_combined_GOterms_sig.*.tsv", full.names = T); sig_files
# note (v6v2)
# the permutation results are still the ones without subsetting the 3 ontogenies
# currently, re-run permuation 
# for now:
# fix this by disgard 'ont' column from the permutation results and 
# use the 'ont' from the custom go annotation files 

out2 = lapply(sig_files, function(s){
  # s = sig_files[6]
  tmp_name = basename(s)
  sp = gsub("(\\w+)_combined_GOterms_sig.*.tsv", "\\1", tmp_name)
  sp = gsub("_1", "", sp)
  print(sp)
  goList <- read_tsv(s, col_names = T) %>%  #, col_types = c('c', 'n', 'c', 'n', 'c')
    mutate(geneID = as.character(geneID)) %>% 
    # dplyr::select(ID) %>% 
    bind_cols('group' = sp)
}) 
out3 = out2 %>% bind_rows() %>%
  dplyr::rename(go_id = ID, ont = subontology) 

#### annotate GO IDs with functional annotation
# way1: get it from the custom go annotation
sel_custom_go_anno = go_annoC %>% 
  dplyr::select(go_id, goTerm, goTerm_def) %>% 
  distinct()
str(out3, 2)
str(sel_custom_go_anno, 2)

go_anno = out3 %>% left_join(sel_custom_go_anno, by = c("go_id" = "go_id")) %>% 
  arrange(group, ont) %>% distinct() %>% 
  group_by(group, ont) %>% 
  mutate(padj = p.adjust(pVal_enrich, "BH")) %>% # 0926
  dplyr::select(group, ont, go_id, goTerm, pVal_enrich, padj, 
                target_in, target_total, bg_in, bg_total, geneID, goTerm_def) # 0926
dim(go_anno) # 466>_2>456   12
go_anno %>% filter(is.na(goTerm)) %>% dplyr::select(go_id) %>% distinct() # nothing is good

# rename sunbirds_flowerpecker
unique(go_anno$group)
go_anno$group = gsub("sunbirds_flowerpecker", "sunbirds", go_anno$group)

# 0926
sum(go_anno$pVal_enrich < 0.05) # 466>_2>456 >> 5968
sum(go_anno$padj < 0.05) # 466>_2>456 >> 5968 # same as not adjusted, so maybe Sara ignore the multiple comparison adjustment...?
# enough permutation will make p value adjust less of a impact
# https://stats.stackexchange.com/questions/93827/permutation-tests-and-multiple-testing
plot(go_anno$padj, go_anno$pVal_enrich)
hist(go_anno$pVal_enrich)
hist(go_anno$padj)
# also, plotting -log10(padj), absolutely weird [GO_eachclade_Hum_2022-09-26_padj.pdf]

path = paste0(out_path, "Final_GO_combined_sig_mart_anno_", Sys.Date(), ".tsv")
# write_tsv(go_anno, path, col_names = T, na = "")
path = paste0(out_path, "Final_GO_combined_sig_mart_anno_2022-11-29.tsv")
path = paste0(out_path, "Final_GO_combined_sig_mart_anno_2023-03-18.tsv") # include _only
path = paste0(out_path, "Final_GO_combined_sig_mart_anno_2023-04-11.tsv") # updated gal6 anno, not include _only
# go_anno = read_tsv(path)

unique(go_anno$group)

x = go_anno %>% 
  filter(group %in% c("honeyeaters", "hummingbirds", "parrots", "sunbirds")) %>% 
  group_by(go_id) %>% 
  mutate(n_clade = length(unique(group)))

go_anno_splitgene = go_anno %>% 
  mutate(lab = paste0(go_id, ":", goTerm), 
         `-log10P` = -log10(pVal_enrich)) %>% 
  mutate(geneIDsplit = strsplit(geneID, split = "/")) %>% 
  unnest(geneIDsplit) %>% 
  dplyr::select(-geneID) %>% dplyr::rename(gene = geneIDsplit)  %>% 
  mutate(lab = fct_reorder(lab, `-log10P`)) %>% 
  distinct()
dim(go_anno_splitgene) # 14690


## load acc and genes ====
### nectar 
## raw acc. cnees (associated with genes (not all acc cnees))
# == a6.1 genes within 100kb (up and downstream) of acc. cnees [in 6_comparelists_20221005.R]
# target_file = list.files(paste0(project_path, "/../linear/"), 
#                          "cnee_ncbigene100kb_acc_intersect_FullTreeExp_\\w+.bed", full.names = T); target_file
target_file = list.files(paste0(project_path, "/../linear/"), 
                         "cnee_ncbigene100kb_acc_intersect_FullTreeExp_\\w+_2.bed", full.names = T); target_file # updated gal6 gnen symbol
acccnees_100kb_withcnees0 = lapply(target_file, function(s){
  # s = target_file[1]
  tmp_name = basename(s)
  tmp_name = gsub("cnee_ncbigene100kb_acc_intersect_FullTreeExp_", "", tmp_name)
  tmp_name = gsub(".bed", "", tmp_name)
  print(paste0("current: ", tmp_name))
  x = read_tsv(s, col_names = F) %>% 
    filter(X6 == "protein_coding") %>%
    dplyr::select(X4, X12, X5, X7, X8) %>%
    dplyr::rename("gene" = "X4", "id" = "X12", "ncbi" = "X5", "accel" = "X7", "total" = "X8") %>% # updated gal6
    # separate(combo, into = c(NA, "ncbi"), sep = "GeneID:", remove = F) %>%  
    bind_cols('group' = tmp_name) %>% 
    dplyr::select(group, gene, ncbi, id, accel, total, id)
}) %>% bind_rows()
head(acccnees_100kb_withcnees0)
# group          gene   ncbi   id        accel total
# <chr>          <chr>  <chr>  <chr>     <dbl> <dbl>
# 1 2wayconvergent KLHL1  418799 CNEE51881     1    44
# 2 2wayconvergent PLXNA4 427941 CNEE501       1   353

length(unique(acccnees_100kb_withcnees0$gene)) # 3523 gene symbols >>5699
length(unique(acccnees_100kb_withcnees0$ncbi)) # 3523 entrez >> 5698

# rm _2 from group
unique(acccnees_100kb_withcnees0$group)
acccnees_100kb_withcnees0 = acccnees_100kb_withcnees0 %>% 
  mutate(group = gsub("_2", "", group))
unique(acccnees_100kb_withcnees0$group)

# rename sunbirds_flowerpecker to sunbirds
unique(acccnees_100kb_withcnees0$group)
acccnees_100kb_withcnees0 = acccnees_100kb_withcnees0 %>% 
  mutate(group = gsub("sunbirds_flowerpecker", "sunbirds", group))
unique(acccnees_100kb_withcnees0$group)

## keep only sets of 4 cloades
acccnees_100kbgg_withcnees_sel = acccnees_100kb_withcnees0 %>% 
  dplyr::rename(clade = group) %>% 
  filter(clade %in% c("honeyeaters", "hummingbirds", "parrots", "sunbirds" ))

acccnees_100kbgg_withcnees_selb = acccnees_100kbgg_withcnees_sel %>% 
  dplyr::select(clade, gene) %>% distinct() 
unique(acccnees_100kbgg_withcnees_selb$clade)

# with clade
go_anno_splitgene_clade = go_anno_splitgene %>% 
  left_join(acccnees_100kbgg_withcnees_selb) %>% 
  filter(!is.na(clade))
unique(go_anno_splitgene_clade$clade)
dim(go_anno_splitgene_clade) # 19458

go_anno_splitgene_clade$clade = factor(go_anno_splitgene_clade$clade, 
                                       levels = c("hummingbirds", "parrots", 
                                                  "honeyeaters", "sunbirds"))
# go_anno_splitgene_clade: look for observed genes from go terms from the 4 clades of cnee-associated genes
# in the case of clade-specific, this will also look for genes from other clades

# compare clade lists ====
clade_ls = lapply(c("hummingbirds", "parrots", "honeyeaters", "sunbirds"),
                  function(s){  tmp = acccnees_100kbgg_withcnees_selb %>% filter(clade == s) %>% pull(gene) })
names(clade_ls) = c("hummingbirds", "parrots", "honeyeaters", "sunbirds")

gene_sim = lapply(c("BP", "CC", "MF"), function(s){
  hsGO <- GOSemSim::godata('org.Hs.eg.db', ont=s, keytype = "SYMBOL")
  GOSemSim::mclusterSim(clade_ls, semData=hsGO, measure="Wang", combine="BMA")
})
# BP
# hummingbirds parrots honeyeaters sunbirds
# hummingbirds                 1.000   0.954       0.919                 0.916
# parrots                      0.954   1.000       0.923                 0.913
# honeyeaters                  0.919   0.923       1.000                 0.904
# sunbirds        0.916   0.913       0.904                 1.000
# 
# CC
# hummingbirds parrots honeyeaters sunbirds
# hummingbirds                 1.000   0.991       0.983                 0.983
# parrots                      0.991   1.000       0.984                 0.984
# honeyeaters                  0.983   0.984       1.000                 0.982
# sunbirds        0.983   0.984       0.982                 1.000
# 
# MF
# hummingbirds parrots honeyeaters sunbirds
# hummingbirds                 1.000   0.968       0.948                 0.951
# parrots                      0.968   1.000       0.949                 0.948
# honeyeaters                  0.948   0.949       1.000                 0.948
# sunbirds        0.951   0.948       0.948                 1.000



### plot ====
## 1. each clade ====
go_anno_tmp = go_anno %>%
  filter(group %in% c("hummingbirds", "parrots", 
                      "honeyeaters", "sunbirds")) %>% 
  mutate(lab = paste0(go_id, ":", goTerm), 
         `-log10P` = -log10(pVal_enrich)) %>% 
  mutate(lab = fct_reorder(lab, `-log10P`))

go_anno_tmp$group = factor(go_anno_tmp$group, levels = c("hummingbirds", "parrots", 
                                                         "honeyeaters", "sunbirds"))
go_anno_tmp$ont = factor(go_anno_tmp$ont, levels = c("biological_process", "molecular_function", "cellular_component"))

pal = brewer.pal(n = 4, name = "Dark2")

pal1 = pal[1]
p1 = ggplot(go_anno_tmp %>% filter(group == "hummingbirds")%>% 
              mutate(lab = fct_reorder(lab, `-log10P`)) ) +
  geom_col(aes(x = `-log10P`, y = lab, fill = group), alpha = 0.7) +
  scale_x_continuous(expand = c(0,0,0.05,0)) +
  facet_grid(facets = group+ont~., space = 'free', scale = 'free') + 
  scale_fill_manual(values = pal1) +
  theme_m + 
  theme(axis.title.y = element_blank(), 
        legend.position = "bottom"); p1

graph_path = paste0(out_path, "GO_eachclade_Hum_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*5, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p1
dev.off()
openfile(graph_path)

pal2 = pal[-1]
p2 = ggplot(go_anno_tmp %>% filter(group != "hummingbirds")  %>% 
              mutate(lab = fct_reorder(lab, `-log10P`)) ) +
  geom_col(aes(x = `-log10P`, y = lab, fill = group), alpha = 0.7) +
  scale_x_continuous(expand = c(0,0,0.05,0)) +
  facet_grid(facets = group+ont~., space = 'free', scale = 'free') + 
  scale_fill_manual(values = pal2) +
  theme_m + 
  theme(axis.title.y = element_blank(), 
        legend.position = "bottom"); p2

graph_path = paste0(out_path, "GO_eachclade_noHum_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*2.75, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p2
dev.off()
openfile(graph_path)

# p = arrangeGrob(grobs = list(p1, p2), nrow = 1)
# grid.draw(p)
# 
# graph_path = paste0(out_path, "GO_eachclade_", Sys.Date(), ".pdf")
# pdf(graph_path, width = Width_HalfCol*4, height = Width_HalfCol*6, pointsize = AxisTxFontSizeSize, onefile = TRUE)
# grid.draw(p)
# dev.off()
# openfile(graph_path)



# convergent ====
go_anno_tmp = go_anno %>%
  filter(group %in% c("hummingbirds", "parrots", 
                      "honeyeaters", "sunbirds")) %>% 
  mutate(lab = paste0(go_id, ":", goTerm), 
         `-log10P` = -log10(pVal_enrich)) %>% 
  mutate(lab = fct_reorder(lab, `-log10P`))

go_anno_tmp$group = factor(go_anno_tmp$group, levels = c("hummingbirds", "parrots", 
                                                         "honeyeaters", "sunbirds"))
go_anno_tmp$ont = factor(go_anno_tmp$ont, levels = c("biological_process", "molecular_function", "cellular_component"))

go_anno_tmp = go_anno_tmp %>% 
  group_by(go_id) %>% 
  mutate(n_clades_byterm = n_distinct(group))

## all sig
go_anno_tmp_conv = go_anno_tmp %>% 
  complete(group, nesting(lab, ont )) %>% 
  # filter(n_clades_byterm >= 2) %>% 
  arrange(n_clades_byterm) 

y_lab = unique(go_anno_tmp_conv$lab)
go_anno_tmp_conv$lab= factor(go_anno_tmp_conv$lab, y_lab)

# size_max = 3

pal = brewer.pal(9, 'Purples')
p1 = ggplot(go_anno_tmp_conv ) +
  geom_tile(aes(x = group, y = lab, fill = `-log10P`), color = 'white', na.rm = FALSE) +
  # scale_x_continuous(expand = c(0,0,0.05,0)) +
  # scale_y_discrete(expand = c(0,0,0,0)) +
  facet_grid(facets = ont~., space = 'free', scale = 'free') +
  # scale_fill_manual(values = pal1) +
  scale_fill_gradient(low = pal[3], high = pal[9], na.value = 'lightyellow') +
  # scale_fill_brewer('seq') +
  labs(x=expression(-~log[10]~P)) +
  theme_m + 
  theme(axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = AxisTxFontSizeSize),
        legend.position = "right"); p1

graph_path = paste0(out_path, "go_4sets_heatmap_all_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2.5, height = Width_HalfCol*7*2, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p1
dev.off()
openfile(graph_path)

# prep comma-sep for katya's python script ====
x = go_anno_tmp_conv %>% group_by(ont) %>% 
  distinct(go_id) %>% 
  mutate(comma = paste0(go_id, collapse = ",")) %>% 
  mutate(space = paste0(go_id, collapse = " ")) %>% 
  dplyr::select(-go_id) %>% 
  distinct()

# x$comma # no need
x$space

# [1] "GO:0000038 GO:0000086 GO:0000122 GO:0000244 GO:0000290 GO:0000301 GO:0000338 GO:0000353 GO:0000395 GO:0000470 GO:0000724 GO:0000731 GO:0001503 GO:0001522 GO:0001541 GO:0001569 GO:0001654 GO:0001755 GO:0001768 GO:0001833 GO:0001886 GO:0001893 GO:0001937 GO:0001944 GO:0001974 GO:0002092 GO:0002244 GO:0002250 GO:0002312 GO:0002323 GO:0002437 GO:0002690 GO:0002934 GO:0003139 GO:0003143 GO:0003360 GO:0005977 GO:0006164 GO:0006282 GO:0006310 GO:0006336 GO:0006357 GO:0006378 GO:0006473 GO:0006508 GO:0006559 GO:0006644 GO:0006686 GO:0006691 GO:0006699 GO:0006730 GO:0006796 GO:0006811 GO:0006812 GO:0006816 GO:0006826 GO:0006851 GO:0006883 GO:0006952 GO:0006954 GO:0006955 GO:0006978 GO:0007005 GO:0007035 GO:0007076 GO:0007157 GO:0007158 GO:0007188 GO:0007193 GO:0007204 GO:0007224 GO:0007252 GO:0007275 GO:0007281 GO:0007338 GO:0007398 GO:0007413 GO:0007417 GO:0007530 GO:0007599 GO:0007605 GO:0007612 GO:0007631 GO:0008206 GO:0008272 GO:0008277 GO:0008283 GO:0008347 GO:0008643 GO:0008652 GO:0009134 GO:0009308 GO:0009416 GO:0009451 GO:0009968 GO:0010466 GO:0010524 GO:0010572 GO:0010587 GO:0010596 GO:0010738 GO:0010745 GO:0010761 GO:0010800 GO:0010801 GO:0010820 GO:0010827 GO:0010838 GO:0010842 GO:0010867 GO:0010882 GO:0010883 GO:0010951 GO:0010952 GO:0014033 GO:0014037 GO:0014824 GO:0015015 GO:0015711 GO:0015807 GO:0015825 GO:0015886 GO:0016024 GO:0016043 GO:0016055 GO:0016125 GO:0016246 GO:0016477 GO:0017121 GO:0017183 GO:0018095 GO:0018146 GO:0019082 GO:0019233 GO:0019370 GO:0019371 GO:0019432 GO:0021559 GO:0021762 GO:0021772 GO:0021819 GO:0021930 GO:0021987 GO:0023052 GO:0030030 GO:0030097 GO:0030099 GO:0030148 GO:0030316 GO:0030322 GO:0030432 GO:0030500 GO:0030502 GO:0030511 GO:0030522 GO:0030705 GO:0030728 GO:0030948 GO:0030953 GO:0031063 GO:0031103 GO:0031113 GO:0031123 GO:0031146 GO:0031167 GO:0031284 GO:0031623 GO:0032060 GO:0032091 GO:0032228 GO:0032308 GO:0032355 GO:0032388 GO:0032482 GO:0032526 GO:0032781 GO:0032868 GO:0032933 GO:0032968 GO:0033077 GO:0033146 GO:0033169 GO:0033173 GO:0033198 GO:0033280 GO:0033504 GO:0033627 GO:0033993 GO:0034405 GO:0034656 GO:0034695 GO:0034773 GO:0034775 GO:0035093 GO:0035239 GO:0035330 GO:0035481 GO:0035493 GO:0035556 GO:0035563 GO:0035590 GO:0035641 GO:0035860 GO:0036035 GO:0038031 GO:0038066 GO:0039694 GO:0040014 GO:0042127 GO:0042270 GO:0042326 GO:0042491 GO:0042593 GO:0042761 GO:0043068 GO:0043113 GO:0043154 GO:0043252 GO:0043279 GO:0043393 GO:0043405 GO:0043457 GO:0043473 GO:0043517 GO:0043518 GO:0043542 GO:0043584 GO:0043586 GO:0043589 GO:0043697 GO:0043923 GO:0043950 GO:0044782 GO:0044790 GO:0045010 GO:0045124 GO:0045444 GO:0045596 GO:0045598 GO:0045651 GO:0045666 GO:0045744 GO:0045762 GO:0045786 GO:0045793 GO:0045820 GO:0045839 GO:0045842 GO:0045860 GO:0045987 GO:0046339 GO:0046548 GO:0046688 GO:0046716 GO:0046839 GO:0048016 GO:0048025 GO:0048148 GO:0048172 GO:0048255 GO:0048341 GO:0048468 GO:0048598 GO:0048643 GO:0048645 GO:0048665 GO:0048675 GO:0048699 GO:0048706 GO:0048712 GO:0048790 GO:0048843 GO:0048845 GO:0048935 GO:0050776 GO:0050807 GO:0050808 GO:0050850 GO:0050909 GO:0051014 GO:0051036 GO:0051123 GO:0051252 GO:0051281 GO:0051298 GO:0051302 GO:0051451 GO:0051452 GO:0051560 GO:0051592 GO:0051593 GO:0051873 GO:0051894 GO:0051918 GO:0051928 GO:0051972 GO:0051974 GO:0055086 GO:0055119 GO:0060070 GO:0060080 GO:0060088 GO:0060100 GO:0060173 GO:0060255 GO:0060261 GO:0060291 GO:0060306 GO:0060315 GO:0060323 GO:0060348 GO:0060349 GO:0060384 GO:0060402 GO:0060411 GO:0060509 GO:0060548 GO:0060612 GO:0060644 GO:0060964 GO:0061003 GO:0061029 GO:0061053 GO:0061077 GO:0061314 GO:0061384 GO:0065008 GO:0070102 GO:0070166 GO:0070168 GO:0070292 GO:0070301 GO:0070306 GO:0070307 GO:0070371 GO:0070373 GO:0070932 GO:0070934 GO:0070986 GO:0071243 GO:0071257 GO:0071318 GO:0071356 GO:0071360 GO:0071361 GO:0071374 GO:0071392 GO:0071455 GO:0071468 GO:0071577 GO:0071625 GO:0071726 GO:0071786 GO:0071901 GO:0071918 GO:0072091 GO:0072711 GO:0080090 GO:0086002 GO:0086064 GO:0086091 GO:0090129 GO:0090136 GO:0090166 GO:0090314 GO:0090666 GO:0097009 GO:0097120 GO:0098703 GO:0098712 GO:0098761 GO:0099175 GO:0120162 GO:0140331 GO:1900181 GO:1900407 GO:1900747 GO:1901137 GO:1901216 GO:1901844 GO:1902166 GO:1902902 GO:1903799 GO:1904045 GO:1904117 GO:1904383 GO:1904659 GO:1904706 GO:1904753 GO:1905665 GO:1905820 GO:1990542 GO:2000036 GO:2000096 GO:2000114 GO:2000310 GO:2000343 GO:2000369 GO:2000727 GO:2000766 GO:2000772 GO:2001020 GO:2001034 GO:2001236 GO:0002931 GO:0006367 GO:0006486 GO:0006897 GO:0007127 GO:0007416 GO:0008219 GO:0010634 GO:0016241 GO:0016266 GO:0019722 GO:0030178 GO:0030336 GO:0031648 GO:0032259 GO:0032967 GO:0043323 GO:0048477 GO:0048813 GO:2000300 GO:0043507"
# [2] "GO:0000138 GO:0000139 GO:0000172 GO:0000220 GO:0000793 GO:0001772 GO:0005665 GO:0005682 GO:0005683 GO:0005685 GO:0005688 GO:0005743 GO:0005770 GO:0005789 GO:0005813 GO:0005838 GO:0005875 GO:0005876 GO:0008290 GO:0014069 GO:0014701 GO:0016323 GO:0016581 GO:0016592 GO:0019005 GO:0022627 GO:0030018 GO:0030133 GO:0030136 GO:0030140 GO:0030175 GO:0030315 GO:0030424 GO:0030425 GO:0030532 GO:0031224 GO:0031253 GO:0031301 GO:0032426 GO:0032432 GO:0032839 GO:0033179 GO:0033391 GO:0034045 GO:0036038 GO:0036064 GO:0036464 GO:0042470 GO:0042564 GO:0042584 GO:0043204 GO:0043679 GO:0044300 GO:0044853 GO:0045178 GO:0046540 GO:0048786 GO:0048788 GO:0055037 GO:0060077 GO:0070847 GO:0071004 GO:0071006 GO:0071007 GO:0071203 GO:0071782 GO:0072562 GO:0089717 GO:0090571 GO:0097228 GO:0097431 GO:0098685 GO:0099060 GO:0110016 GO:1904115 GO:1990909 GO:0032580 GO:0098793 GO:0098794"        
# [3] "GO:0000900 GO:0001094 GO:0001537 GO:0001540 GO:0001614 GO:0003688 GO:0003697 GO:0003713 GO:0004065 GO:0004089 GO:0004222 GO:0004252 GO:0004364 GO:0004439 GO:0004652 GO:0004712 GO:0004867 GO:0004888 GO:0004931 GO:0004972 GO:0005044 GO:0005096 GO:0005114 GO:0005123 GO:0005154 GO:0005184 GO:0005385 GO:0005522 GO:0005537 GO:0005545 GO:0005546 GO:0005547 GO:0008035 GO:0008194 GO:0008195 GO:0008233 GO:0008237 GO:0008301 GO:0008494 GO:0008499 GO:0008526 GO:0015075 GO:0015194 GO:0015232 GO:0015278 GO:0015279 GO:0015347 GO:0016004 GO:0016018 GO:0016208 GO:0016493 GO:0016504 GO:0016594 GO:0016709 GO:0016740 GO:0016836 GO:0017111 GO:0017128 GO:0017147 GO:0019213 GO:0019828 GO:0019838 GO:0019855 GO:0019887 GO:0019957 GO:0022841 GO:0022857 GO:0023026 GO:0030169 GO:0030331 GO:0030374 GO:0030414 GO:0031072 GO:0031434 GO:0031681 GO:0032454 GO:0035005 GO:0035091 GO:0036002 GO:0042277 GO:0043022 GO:0043225 GO:0043325 GO:0043621 GO:0045504 GO:0046790 GO:0046873 GO:0047617 GO:0048306 GO:0050321 GO:0050750 GO:0050839 GO:0051015 GO:0051020 GO:0051059 GO:0051378 GO:0051721 GO:0051787 GO:0052689 GO:0061575 GO:0070006 GO:0070063 GO:0070181 GO:0072542 GO:0090556 GO:0097677 GO:0097718 GO:0102991 GO:0140439 GO:0140658 GO:1990825 GO:0002020 GO:0004497 GO:0008168"    

### converg
go_anno_tmp_conv = go_anno_tmp %>% 
  filter(n_clades_byterm >= 2) %>% 
  arrange(n_clades_byterm, `-log10P`) %>% 
  complete(group, nesting(lab ))
  
y_lab = unique(go_anno_tmp_conv$lab)
go_anno_tmp_conv$lab= factor(go_anno_tmp_conv$lab, y_lab)

# size_max = 3

pal = brewer.pal(9, 'Purples')
p1 = ggplot(go_anno_tmp_conv ) +
  geom_tile(aes(x = group, y = lab, fill = `-log10P`), color = 'white', na.rm = FALSE) +
  # scale_x_continuous(expand = c(0,0,0.05,0)) +
  # scale_y_discrete(expand = c(0,0,0,0)) +
  # facet_grid(facets = ont~., space = 'free', scale = 'free') +
  # scale_fill_manual(values = pal1) +
  scale_fill_gradient(low = pal[3], high = pal[9], na.value = 'lightyellow') +
  # scale_fill_brewer('seq') +
  labs(x=expression(-~log[10]~P)) +
  theme_m + 
  theme(axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = AxisTxFontSizeSize),
        legend.position = "right"); p1

graph_path = paste0(out_path, "go_4sets_heatmap_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*1.25, height = Width_HalfCol*1.5, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p1
dev.off()
openfile(graph_path)



# 2. 4sets only ====
## cnee convergent 
unique(go_anno$group)

go_anno_tmp = go_anno %>%
  filter(group %in% c("4sets")) %>% 
  mutate(lab = paste0(go_id, ":", goTerm), 
         `-log10P` = -log10(pVal_enrich)) %>% 
  mutate(lab = fct_reorder(lab, `-log10P`))
unique(go_anno_tmp$group)

go_anno_tmp$group = factor(go_anno_tmp$group, levels = c("4sets"))
go_anno_tmp$ont = factor(go_anno_tmp$ont, levels = c("biological_process", "molecular_function", "cellular_component"))

path = paste0(out_path, "4sets_sigGO_", Sys.Date(), ".tsv")
# write_tsv(go_anno_tmp, path, col_names = T, na = "")


# pal = brewer.pal(n = 4, name = "Dark2")
# pal1 = pal[1]

pal = brewer.pal(n = 11, name = "Spectral")
pal1 = pal[2]
p1 = ggplot(go_anno_tmp ) +
  geom_col(aes(x = `-log10P`, y = lab, fill = group), alpha = 0.7) +
  scale_x_continuous(expand = c(0,0,0.05,0)) +
  labs(x=expression(-~log[10]~P)) +
  facet_grid(facets = group+ont~., space = 'free', scale = 'free') +
  scale_fill_manual(values = pal1) +
  theme_m + 
  theme(axis.title.y = element_blank(), 
        legend.position = "bottom"); p1

graph_path = paste0(out_path, "GO_4sets_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2.3, height = Width_HalfCol*4, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p1
dev.off()
openfile(graph_path)

# heatmap
pal = brewer.pal(9, 'Reds')
p1 = ggplot(go_anno_tmp ) +
  geom_tile(aes(x = 1, y = lab, fill = `-log10P`), color = 'white') +
  scale_x_continuous(expand = c(0,0,0.05,0)) +
  scale_y_discrete(expand = c(0,0,0,0)) +
  facet_grid(facets = ont~., space = 'free', scale = 'free') +
  # scale_fill_manual(values = pal1) +
  scale_fill_gradient(low = pal[2], high = "#dc2943ff") +
  labs(x=expression(-~log[10]~P)) +
  theme_m + 
  theme(axis.title.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = AxisTxFontSizeSize),
        legend.position = "right"); p1

graph_path = paste0(out_path, "go_4sets_heatmap_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*3, height = Width_HalfCol*4.5, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p1
dev.off()
openfile(graph_path)

# 4sets BP
# pal1 = pal[1]
p1 = ggplot(go_anno_tmp %>% # filter(group == "4sets") %>% 
              filter(ont == "biological_process") %>% 
              mutate(lab = fct_reorder(lab, `-log10P`)) ) +
  geom_col(aes(x = `-log10P`, y = lab, fill = group), alpha = 0.7) +
  scale_x_continuous(expand = c(0,0,0.05,0)) +
  labs(x=expression(-~log[10]~P)) +
  # facet_grid(facets = group+ont~., space = 'free', scale = 'free') +
  scale_fill_manual(values = pal1) +
  theme_m + 
  theme(axis.title.y = element_blank(), 
        legend.position = "bottom"); p1

graph_path = paste0(out_path, "GO_4sets_bp_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*1.75, height = Width_HalfCol*3, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p1
dev.off()
openfile(graph_path)

# 4sets BP heatmap
pal = brewer.pal(n = 11, name = "Spectral")
pal1 = pal[2]
tmp_4sets_bp = go_anno_tmp %>% filter(group == "4sets") %>% filter(ont == "biological_process") %>% 
  mutate(lab = fct_reorder(lab, `-log10P`))

path = paste0(out_path, "4sets_bp_sigGO_", Sys.Date(), ".tsv")
# write_tsv(tmp_4sets_bp, path, col_names = T, na = "")

pal = brewer.pal(9, 'Reds')
p1 = ggplot( tmp_4sets_bp ) +
  geom_tile(aes(x = 1, y = lab, fill = `-log10P`), color = 'white') +
  scale_x_continuous(expand = c(0,0,0.05,0)) +
  scale_y_discrete(expand = c(0,0,0,0)) +
  # facet_grid(facets = group+ont~., space = 'free', scale = 'free') +
  # scale_fill_manual(values = pal1) +
  scale_fill_gradient(low = pal[2], high = "#dc2943ff") +
  labs(x=expression(-~log[10]~P)) +
  theme_m + 
  theme(axis.title.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = AxisTxFontSizeSize),
        legend.position = "right"); p1

graph_path = paste0(out_path, "GO_4sets_bp_heatmap_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*1.75, height = Width_HalfCol*3, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p1
dev.off()
openfile(graph_path)


# 4sets BP top20 -log10P ====
tmp_4sets_bp_top20 = go_anno_tmp %>% filter(group == "4sets") %>% 
  filter(ont == "biological_process") %>% 
  slice_max(`-log10P`, n = 20) %>% 
  mutate(lab = fct_reorder(lab, `-log10P`))


# heatmap
pal = brewer.pal(9, 'Reds')
p1 = ggplot( tmp_4sets_bp_top20 ) +
  geom_tile(aes(x = 1, y = lab, fill = `-log10P`), color = 'white') +
  scale_x_continuous(expand = c(0,0,0.05,0)) +
  scale_y_discrete(expand = c(0,0,0,0)) +
  # facet_grid(facets = group+ont~., space = 'free', scale = 'free') +
  # scale_fill_manual(values = pal1) +
  scale_fill_gradient(low = pal[2], high = "#dc2943ff") +
  labs(x=expression(-~log[10]~P)) +
  theme_m + 
  theme(axis.title.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = AxisTxFontSizeSize),
        legend.position = "right"); p1

graph_path = paste0(out_path, "GO_4sets_bp_top20P_heatmap_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*1.425, height = Width_HalfCol*0.8, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p1
dev.off()
openfile(graph_path)

tmp_4sets_bp_top20_splitgene = tmp_4sets_bp_top20 %>% 
  mutate(geneIDsplit = strsplit(geneID, split = "/")) %>% 
  unnest(geneIDsplit) %>% 
  dplyr::select(-geneID) %>% dplyr::rename(gene = geneIDsplit)  %>% 
  mutate(lab = fct_reorder(lab, `-log10P`))

# p2 = ggplot( tmp_4sets_bp_top20_splitgene ) +
#   geom_tile(aes(x = gene, y = lab), color = 'white', fill = pal[5]) +
#   # scale_x_continuous(expand = c(0,0,0.05,0)) +
#   # scale_y_discrete(expand = c(0,0,0,0)) +
#   # facet_grid(facets = group+ont~., space = 'free', scale = 'free') +
#   # scale_fill_manual(values = pal1) +
#   # scale_fill_gradient(low = pal[2], high = pal[9]) +
#   labs(x=expression(-~log[10]~P)) +
#   theme_m + 
#   theme(axis.title.y = element_blank(), 
#         # axis.text.x = element_blank(),
#         # axis.ticks.x = element_blank(),
#         legend.position = "none"); p2

# 4sets BP top10 -log10P ====
tmp_4sets_bp_top10 = go_anno_tmp %>% filter(group == "4sets") %>% 
  filter(ont == "biological_process") %>% 
  slice_max(`-log10P`, n = 10) %>% 
  mutate(lab = fct_reorder(lab, `-log10P`))

# heatmap
pal = brewer.pal(9, 'Reds')
p1 = ggplot( tmp_4sets_bp_top10 ) +
  geom_tile(aes(x = 1, y = lab, fill = `-log10P`), color = 'white') +
  scale_x_continuous(expand = c(0,0,0.05,0)) +
  scale_y_discrete(expand = c(0,0,0,0)) +
  # facet_grid(facets = group+ont~., space = 'free', scale = 'free') +
  # scale_fill_manual(values = pal1) +
  scale_fill_gradient(low = pal[2], high = "#dc2943ff") +
  labs(x=expression(-~log[10]~P)) +
  theme_m + 
  theme(axis.title.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = AxisTxFontSizeSize),
        legend.position = "right"); p1

graph_path = paste0(out_path, "GO_4sets_bp_top10P_heatmap_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*1.1, height = Width_HalfCol*0.5, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p1
dev.off()
openfile(graph_path)

# 4sets BP top30 -log10P ====
tmp_4sets_bp_top30 = go_anno_tmp %>% filter(group == "4sets") %>% 
  filter(ont == "biological_process") %>% 
  slice_max(`-log10P`, n = 30) %>% 
  mutate(lab = fct_reorder(lab, `-log10P`))

# heatmap
pal = brewer.pal(9, 'Reds')
p1 = ggplot( tmp_4sets_bp_top30 ) +
  geom_tile(aes(x = 1, y = lab, fill = `-log10P`), color = 'white') +
  scale_x_continuous(expand = c(0,0,0.05,0)) +
  scale_y_discrete(expand = c(0,0,0,0)) +
  # facet_grid(facets = group+ont~., space = 'free', scale = 'free') +
  # scale_fill_manual(values = pal1) +
  scale_fill_gradient(low = pal[2], high = "#dc2943ff") +
  labs(x=expression(-~log[10]~P)) +
  theme_m + 
  theme(axis.title.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = AxisTxFontSizeSize),
        legend.position = "right"); p1

graph_path = paste0(out_path, "GO_4sets_bp_top30P_heatmap_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*1.425, height = Width_HalfCol*0.9, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p1
dev.off()
openfile(graph_path)

# 4sets BP top20 termsim ====
tmp_4sets_bp_top20term = go_anno_tmp %>% filter(group == "4sets") %>% 
  filter(go_id %in% top20_termsim) %>% 
  mutate(lab = goTerm) %>% 
  mutate(lab = fct_reorder(lab, `-log10P`))

pal = brewer.pal(9, 'Reds')
p1 = ggplot( tmp_4sets_bp_top20term ) +
  geom_tile(aes(x = 1, y = lab, fill = `-log10P`), color = 'white') +
  scale_x_continuous(expand = c(0,0,0.05,0)) +
  scale_y_discrete(expand = c(0,0,0,0)) +
  # facet_grid(facets = group+ont~., space = 'free', scale = 'free') +
  # scale_fill_manual(values = pal1) +
  scale_fill_gradient(low = pal[2], high = "#dc2943ff") +
  labs(x=expression(-~log[10]~P)) +
  theme_m + 
  theme(axis.title.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = AxisTxFontSizeSize),
        legend.position = "right"); p1

graph_path = paste0(out_path, "GO_4sets_bp_top20termsim_heatmap_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*1.2, height = Width_HalfCol*0.8, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p1
dev.off()
openfile(graph_path)


# 4sets BP top10 termsim ====
tmp_4sets_bp_top10term = go_anno_tmp %>% filter(group == "4sets") %>% 
  filter(go_id %in% top10_termsim) %>% 
  mutate(lab = goTerm) %>% 
  mutate(lab = fct_reorder(lab, `-log10P`))

pal = brewer.pal(9, 'Reds')
p1 = ggplot( tmp_4sets_bp_top10term ) +
  geom_tile(aes(x = 1, y = lab, fill = `-log10P`), color = 'white') +
  scale_x_continuous(expand = c(0,0,0.05,0)) +
  scale_y_discrete(expand = c(0,0,0,0)) +
  # facet_grid(facets = group+ont~., space = 'free', scale = 'free') +
  # scale_fill_manual(values = pal1) +
  scale_fill_gradient(low = pal[2], high = "#dc2943ff") +
  labs(x=expression(-~log[10]~P)) +
  theme_m + 
  theme(axis.title.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = AxisTxFontSizeSize),
        legend.position = "right"); p1

graph_path = paste0(out_path, "GO_4sets_bp_top10termsim_heatmap_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*0.9, height = Width_HalfCol*0.5, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p1
dev.off()
openfile(graph_path)

# 4sets BP top30 termsim ====
tmp_4sets_bp_top30term = go_anno_tmp %>% filter(group == "4sets") %>% 
  filter(go_id %in% top30_termsim) %>% 
  mutate(lab = goTerm) %>% 
  mutate(lab = fct_reorder(lab, `-log10P`))

pal = brewer.pal(9, 'Reds')
p1 = ggplot( tmp_4sets_bp_top30term ) +
  geom_tile(aes(x = 1, y = lab, fill = `-log10P`), color = 'white') +
  scale_x_continuous(expand = c(0,0,0.05,0)) +
  scale_y_discrete(expand = c(0,0,0,0)) +
  # facet_grid(facets = group+ont~., space = 'free', scale = 'free') +
  # scale_fill_manual(values = pal1) +
  scale_fill_gradient(low = pal[2], high = "#dc2943ff") +
  labs(x=expression(-~log[10]~P)) +
  theme_m + 
  theme(axis.title.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = AxisTxFontSizeSize),
        legend.position = "right"); p1

graph_path = paste0(out_path, "GO_4sets_bp_top30termsim_heatmap_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*1.2, height = Width_HalfCol*0.7, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p1
dev.off()
openfile(graph_path)

#### visual and sound ====
# GO:0001654:eye development 
# GO:0007601:visual perception
# GO:0071456:cellular response to hypoxia
# GO:0048596:embryonic camera type eye morphogenesis
# GO:0007605:sensory perception of sound

x = go_anno_tmp %>% filter(group == "4sets") %>% 
  filter(go_id %in% c("GO:0001654", "GO:0007601", "GO:0071456", "GO:0048596"))


go_anno_splitgene = go_anno %>% 
  mutate(lab = paste0(go_id, ":", goTerm), 
         `-log10P` = -log10(pVal_enrich)) %>% 
  mutate(geneIDsplit = strsplit(geneID, split = "/")) %>% 
  unnest(geneIDsplit) %>% 
  dplyr::select(-geneID) %>% dplyr::rename(gene = geneIDsplit)  %>% 
  mutate(lab = fct_reorder(lab, `-log10P`))

x = go_anno_splitgene %>% filter(group == "4sets") %>% 
  filter(go_id %in% c("GO:0001654", "GO:0007601", "GO:0071456", "GO:0048596", "GO:0007605"))

# with genes ====
go_anno_splitgene_clade_4sets = go_anno_splitgene_clade %>% 
  filter(group == '4sets') %>% 
  mutate(lab = fct_drop(lab))
unique(go_anno_splitgene_clade_4sets$clade)

# v2: num cnee per term per clade
geneCount = go_anno_splitgene_clade_4sets %>% 
  group_by(go_id, clade) %>% # per term per clade, how many (unique) cnees, (sum of all genes)
  mutate(n_genes_perclade_perterm = length(unique(gene))) %>% 
  group_by(go_id) %>% 
  mutate(n_clades = length(unique(clade)),
         sum_gene = sum(n_genes_perclade_perterm)) 

geneCount_b = geneCount %>% 
  dplyr::select(-gene) %>% distinct() %>% 
  group_by(go_id) %>% 
  mutate(n_clades = length(unique(clade)),
         sum_gene = sum(n_genes_perclade_perterm)) 
dim(geneCount_b) # 294  17

path = paste0(out_path, "4sets_sigGO_withgene_", Sys.Date(), ".tsv")
# write_tsv(geneCount_b, path, col_names = T, na = "")
path = paste0(out_path, "4sets_sigGO_withgene_2022-12-12.tsv")
path = paste0(out_path, "4sets_sigGO_withgene_2023-04-11.tsv")

colp = brewer.pal(4, 'Dark2')
# colp2 = c(rep(colp[1], 2), rep(colp[2], 2), colp[3])

N = 3
size_max = ceiling(max(geneCount$n_genes_perclade_perterm)/N)*N
# x = go_anno_gene_sel_cn2 %>% filter(goID =="GO:0051179")

theme_a = theme(axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1),
                axis.text.y = element_text(hjust = 1),
                legend.title = element_text(size = AxisTxFontSizeSize),
                legend.position = "right")

# color of dots in grey ====
# 2.2.1 sort by n_clades, sum_gene ====
y = geneCount_b %>% arrange(n_clades, sum_gene)
y_lab = unique(y$lab)
y$lab = factor(y$lab, y_lab)
y_4sets = y

size_max = 3

# 3
p_term3 = y %>% ggplot() +
  geom_point(aes(y = lab, #fct_reorder2(lab, n_genes_perclade_perterm, sum_gene, .desc = F), 
                 x = clade, size = n_genes_perclade_perterm), color = 'grey50') + # color = clade, 
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_x_discrete(name = "# of CNEEs in accelerated state") +
  scale_color_manual(values = colp, guide = "none") +
  scale_size_continuous(name = "Number of genes", range = c(0.2, size_max)) +
  # guides() +
  theme_m + 
  theme_a; p_term3

# use p_term23's order
# x3 = fct_reorder2(geneCount$lab, geneCount$n_genes_perclade_perterm, geneCount$sum_gene, .desc = F)
# x3 = geneCount %>% group_by(ont) %>% mutate(lab = fct_reorder2(lab, n_genes_perclade_perterm, sum_gene, .desc = F))
# x3_2 = x3$lab

pal = brewer.pal(9, 'Reds')
p3 = ggplot( y ) +
  geom_tile(aes(x = 1, 
                y = lab, #fct_relevel(lab, levels(x3_2)), 
                fill = `-log10P`), color = 'white') +
  scale_x_continuous(expand = c(0,0,0.05,0), breaks = 1, labels = expression(-~log[10]~P)) +
  scale_y_discrete(expand = c(0,0,0,0)) +
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_fill_manual(values = pal1) +
  scale_fill_gradient(low = pal[2], high = "#dc2943ff", name = expression(-~log[10]~P)) +
  # labs(x=expression(-~log[10]~P)) +
  theme_m + 
  theme_a; p3

p = cowplot::plot_grid(p3, p_term3, ncol= 2, scale = rep(c(0.9566691, 1), 3), align = 'v'); 
# rel_widths=c(1.3, 1.3, 1, 1.3)) # labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2)

graph_path = paste0(out_path, "perterm_Ngeneperclade_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2.6*2, height = Width_HalfCol*4.5, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p
dev.off()
openfile(graph_path)

# STable ====
geneCount_stable = geneCount %>% dplyr::select(-group, -padj, -pVal_enrich, -lab) %>% 
  dplyr::rename(ontology = ont) %>% 
  arrange(ontology, desc(n_clades), desc(sum_gene)) %>% 
  dplyr::select(-sum_gene) %>% 
  group_by(go_id, clade) %>% 
  mutate(gene = paste0(gene, collapse = "|")) %>% 
  distinct() %>% 
  dplyr::select(ontology, go_id, goTerm, goTerm_def, 
                target_in:bg_total, `-log10P`:n_clades)

path = paste0(out_path, "4sets_sigGO_withgene_STable_", Sys.Date(), ".tsv")
# write_tsv(geneCount_stable, path, col_names = T, na = "")
path = paste0(out_path, "4sets_sigGO_withgene_STable_2023-01-13.tsv")
path = paste0(out_path, "4sets_sigGO_withgene_STable_2023-04-11.tsv")
# geneCount_stable = read_tsv(path)

geneCount_stable %>% group_by(ontology, n_clades) %>% tally() %>% 
  mutate(sum = sum(n), percent = n/sum)
# ontology           n_clades     n   sum percent
# <chr>                 <dbl> <int> <int>   <dbl>
# 1 biological_process        1    13   186  0.0699
# 2 biological_process        2    30   186  0.161 
# 3 biological_process        3    51   186  0.274 
# 4 biological_process        4    92   186  0.495 
# 5 cellular_component        1     4    60  0.0667
# 6 cellular_component        2     6    60  0.1   
# 7 cellular_component        3     6    60  0.1   
# 8 cellular_component        4    44    60  0.733 
# 9 molecular_function        1     5    48  0.104 
# 10 molecular_function        2     2    48  0.0417
# 11 molecular_function        3     9    48  0.188 
# 12 molecular_function        4    32    48  0.667 

# prep comma-sep for katya's python script ====
x = geneCount %>% group_by(ont) %>% 
  distinct(go_id) %>% 
  mutate(comma = paste0(go_id, collapse = ",")) %>% 
  mutate(space = paste0(go_id, collapse = " ")) %>% 
  dplyr::select(-go_id) %>% 
  distinct()

# x$comma # no need
x$space
# BP
# "GO:0000338 GO:0001516 GO:0001522 GO:0001654 GO:0001764 GO:0001833 GO:0001886 GO:0001893 GO:0001937 GO:0001944 GO:0002011 GO:0002098 GO:0002931 GO:0003376 GO:0005977 GO:0006357 GO:0006367 GO:0006417 GO:0006486 GO:0006605 GO:0006699 GO:0006730 GO:0006816 GO:0006829 GO:0006836 GO:0006851 GO:0006882 GO:0006897 GO:0007158 GO:0007194 GO:0007204 GO:0007252 GO:0007405 GO:0007413 GO:0007548 GO:0007596 GO:0007599 GO:0007605 GO:0007631 GO:0008015 GO:0008033 GO:0008209 GO:0008219 GO:0008652 GO:0009101 GO:0009312 GO:0009451 GO:0010001 GO:0010467 GO:0010596 GO:0010738 GO:0010800 GO:0010820 GO:0010880 GO:0014033 GO:0014912 GO:0015711 GO:0015721 GO:0015804 GO:0015807 GO:0015813 GO:0016043 GO:0016055 GO:0016125 GO:0016266 GO:0016572 GO:0017183 GO:0018146 GO:0019082 GO:0021762 GO:0021819 GO:0030041 GO:0030111 GO:0030178 GO:0030278 GO:0030336 GO:0030500 GO:0030511 GO:0030705 GO:0030728 GO:0031175 GO:0031623 GO:0031648 GO:0032355 GO:0032526 GO:0032786 GO:0032968 GO:0033137 GO:0033147 GO:0033173 GO:0033280 GO:0033627 GO:0035988 GO:0042552 GO:0042754 GO:0042789 GO:0043473 GO:0043507 GO:0043517 GO:0043697 GO:0043966 GO:0045475 GO:0045666 GO:0045793 GO:0048568 GO:0048665 GO:0048706 GO:0048745 GO:0048935 GO:0050905 GO:0051013 GO:0051123 GO:0051560 GO:0051591 GO:0051592 GO:0051918 GO:0060261 GO:0060416 GO:0060840 GO:0061003 GO:0065003 GO:0065008 GO:0070102 GO:0070509 GO:0070542 GO:0070588 GO:0071243 GO:0071320 GO:0071377 GO:0071577 GO:0090043 GO:0090129 GO:0097503 GO:0098703 GO:1900182 GO:1900407 GO:1900747 GO:1901137 GO:1901224 GO:1901387 GO:1902459 GO:1903799 GO:1904262 GO:2000369 GO:2001235"
# CC
# "GO:0000178 GO:0005875 GO:0005901 GO:0008180 GO:0008305 GO:0014069 GO:0014701 GO:0016592 GO:0030139 GO:0030425 GO:0030863 GO:0032420 GO:0032426 GO:0032580 GO:0032584 GO:0042584 GO:0044853 GO:0048471 GO:0048786 GO:0060076 GO:0070847 GO:0072562 GO:0097386 GO:0097440 GO:0098793 GO:0098831 GO:0140672"      
# MF
# "GO:0004181 GO:0004866 GO:0005044 GO:0005123 GO:0005154 GO:0005246 GO:0005262 GO:0005385 GO:0005546 GO:0008195 GO:0008270 GO:0008301 GO:0008324 GO:0008373 GO:0008514 GO:0015179 GO:0015269 GO:0015279 GO:0016407 GO:0016493 GO:0016709 GO:0016779 GO:0016836 GO:0016922 GO:0017147 GO:0019887 GO:0019957 GO:0030159 GO:0030331 GO:0032266 GO:0034185 GO:0034237 GO:0035091 GO:0043325 GO:0046332 GO:0050998 GO:0051015 GO:0051059 GO:0061575 GO:0070006 GO:0070181 GO:0097573"  

# 2.2.2 sort by p values ====
size_max = 3

y = geneCount_b %>%
  arrange(`-log10P`)
# mutate(lab = fct_reorder(lab, `-log10P`)) # doesn't work perfectly
y_lab = unique(y$lab)
y$lab = factor(y$lab, y_lab)

# 3
p_term4 = y %>% ggplot() +
  geom_point(aes(#y = fct_reorder2(lab, n_genes_perclade_perterm, sum_gene, .desc = F), 
                 y = lab,
                 x = clade, size = n_genes_perclade_perterm), color = 'grey50') + # color = clade, 
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_x_discrete(name = "# of CNEEs in accelerated state") +
  scale_color_manual(values = colp, guide = "none") +
  scale_size_continuous(name = "Number of genes", range = c(0.2, size_max)) +
  # guides() +
  theme_m + 
  theme_a; p_term4

# use p_term23's order
# x3 = fct_reorder2(geneCount$lab, geneCount$n_genes_perclade_perterm, geneCount$sum_gene, .desc = F)
# x3 = geneCount %>% group_by(ont) %>% mutate(lab = fct_reorder2(lab, n_genes_perclade_perterm, sum_gene, .desc = F))
# x3_2 = x3$lab

pal = brewer.pal(9, 'Reds')
p4 = ggplot( y ) +
  geom_tile(aes(x = 1, #y = fct_relevel(lab, levels(x3_2)), 
                y = lab, fill = `-log10P`), color = 'white') +
  scale_x_continuous(expand = c(0,0,0.05,0), breaks = 1, labels = expression(-~log[10]~P)) +
  scale_y_discrete(expand = c(0,0,0,0)) +
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_fill_manual(values = pal1) +
  scale_fill_gradient(low = pal[2], high = "#dc2943ff", name = expression(-~log[10]~P)) +
  # labs(x=expression(-~log[10]~P)) +
  theme_m + 
  theme_a; p4

p = cowplot::plot_grid(p4, p_term4, ncol= 2, scale = rep(c(0.9566691, 1), 3), align = 'v'); 
# rel_widths=c(1.3, 1.3, 1, 1.3)) # labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2)

graph_path = paste0(out_path, "perterm_Ngeneperclade_bypval_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*1.65*2, height = Width_HalfCol*2.7, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p
dev.off()
openfile(graph_path)


## plot only 4-way convergent pathways====
size_max = 3

# 3
colp = brewer.pal(4, 'Dark2')
tmp = geneCount_b %>% filter(n_clades >=4) %>% arrange(n_clades, sum_gene)
y_lab = unique(tmp$lab)
tmp$lab = factor(tmp$lab, y_lab)

p_term3 = tmp %>% ggplot() +
  geom_point(aes(y = lab, #fct_reorder2(lab, n_genes_perclade_perterm, sum_gene), 
                 x = clade, size = n_genes_perclade_perterm), color = 'grey50') + # color = clade, 
  # facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_x_discrete(name = "# of CNEEs in accelerated state") +
  scale_color_manual(values = colp, guide = "none") +
  scale_size_continuous(name = "Number of genes", range = c(0.2, size_max)) +
  # guides() +
  theme_m + 
  theme_a; p_term3

# use p_term23's order
# x3 = fct_reorder2(tmp$lab, tmp$n_genes_perclade_perterm, tmp$sum_gene)
pal = brewer.pal(9, 'Reds')
p3 = ggplot(tmp ) +
  geom_tile(aes(x = 1, Ey = fct_relevel(lab, levels(x3)), 
                y = lab,
                fill = `-log10P`), color = 'white') +
  scale_x_continuous(expand = c(0,0,0.05,0), breaks = 1, labels = expression(-~log[10]~P)) +
  scale_y_discrete(expand = c(0,0,0,0)) +
  # facet_grid(facets = group+ont~., space = 'free', scale = 'free') +
  # scale_fill_manual(values = pal1) +
  scale_fill_gradient(low = pal[2], high = "#dc2943ff") +
  # labs(x=expression(-~log[10]~P)) +
  theme_m + 
  theme_a; p3

p = cowplot::plot_grid(p3, p_term3, ncol= 2, scale = rep(c(0.9566691, 1), 1), align = 'v'); # good scale
# rel_widths=c(1.3, 1.3, 1, 1.3)) # labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2)

graph_path = paste0(out_path, "perterm_Ngeneperclade_4_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2*1.6, height = Width_HalfCol*2, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p
dev.off()
openfile(graph_path)


###
pal2 = pal[-1]
p2 = ggplot(go_anno_tmp %>% filter(group != "4sets")%>% 
              mutate(lab = fct_reorder(lab, `-log10P`)) ) +
  geom_col(aes(x = `-log10P`, y = lab, fill = group), alpha = 0.7) +
  scale_x_continuous(expand = c(0,0,0.05,0)) +
  facet_grid(facets = group+ont~., space = 'free', scale = 'free') +
  scale_fill_manual(values = pal2) +
  theme_m + 
  theme(axis.title.y = element_blank(), 
        legend.position = "bottom"); p2

graph_path = paste0(out_path, "GO_2and3way_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2.5, height = Width_HalfCol*2.5, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p2
dev.off()
openfile(graph_path)

## gene-lv convergent 
unique(go_anno$group)

go_anno_tmp = go_anno %>%
  filter(group %in% c("genelv_conv_2way", "genelv_conv_3way", "genelv_conv_4way")) %>% 
  mutate(lab = paste0(go_id, ":", goTerm), 
         `-log10P` = -log10(pVal_enrich)) %>% 
  mutate(lab = fct_reorder(lab, `-log10P`))
unique(go_anno_tmp$group)

go_anno_tmp$group = factor(go_anno_tmp$group, levels = c("genelv_conv_2way", "genelv_conv_3way", "genelv_conv_4way"))
go_anno_tmp$ont = factor(go_anno_tmp$ont, levels = c("biological_process", "molecular_function", "cellular_component"))

pal = brewer.pal(n = 4, name = "Dark2")

pal1 = pal[1]
p1 = ggplot(go_anno_tmp %>% filter(group == "genelv_conv_2way") %>% 
              mutate(lab = fct_reorder(lab, `-log10P`)) ) +
  geom_col(aes(x = `-log10P`, y = lab, fill = group), alpha = 0.7) +
  scale_x_continuous(expand = c(0,0,0.05,0)) +
  facet_grid(facets = group+ont~., space = 'free', scale = 'free') +
  scale_fill_manual(values = pal1) +
  theme_m + 
  theme(axis.title.y = element_blank(), 
        legend.position = "bottom"); p1

graph_path = paste0(out_path, "GO_genelv_conv_2way_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*2.2, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p1
dev.off()
openfile(graph_path)


pal2 = pal[-1]
p2 = ggplot(go_anno_tmp %>% filter( group %in% c("genelv_conv_3way", "genelv_conv_4way")  ) %>% 
              mutate(lab = fct_reorder(lab, `-log10P`)) ) +
  geom_col(aes(x = `-log10P`, y = lab, fill = group), alpha = 0.7) +
  scale_x_continuous(expand = c(0,0,0.05,0)) +
  facet_grid(facets = group+ont~., space = 'free', scale = 'free') +
  scale_fill_manual(values = pal2) +
  theme_m + 
  theme(axis.title.y = element_blank(), 
        legend.position = "bottom"); p2

graph_path = paste0(out_path, "GO_genelv_conv_3and4way_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*1.5, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p2
dev.off()
openfile(graph_path)

# p = arrangeGrob(grobs = list(p1, p2), nrow = 2)
# grid.draw(p)
# 
# graph_path = paste0(out_path, "GO_4sets_2way_", Sys.Date(), ".pdf")
# pdf(graph_path, width = Width_HalfCol*4, height = Width_HalfCol*1.75, pointsize = AxisTxFontSizeSize, onefile = TRUE)
# grid.draw(p)
# dev.off()
# openfile(graph_path)

# 3. clade specific [with genes] ====
go_anno_splitgene_clade_cladesp = go_anno_splitgene_clade %>% 
  filter( grepl("_only", group)) %>% 
  mutate(lab = fct_drop(lab))
unique(go_anno_splitgene_clade_cladesp$clade)

# v2: num cnee per term per clade
geneCount = go_anno_splitgene_clade_cladesp %>% 
  group_by(go_id, clade) %>% # per term per clade, how many (unique) cnees, (sum of all genes)
  mutate(n_genes_perclade_perterm = length(unique(gene))) %>% 
  group_by(go_id) %>% 
  mutate(n_clades = length(unique(clade)),
         sum_gene = sum(n_genes_perclade_perterm)) 

geneCount_b = geneCount %>% 
  dplyr::select(-gene) %>% distinct() %>% 
  group_by(go_id) %>% 
  mutate(n_clades = length(unique(clade)),
         sum_gene = sum(n_genes_perclade_perterm)) 
dim(geneCount_b) # 676  17

path = paste0(out_path, "4sets_sigGO_withgene_cladespecific_", Sys.Date(), ".tsv")
# write_tsv(geneCount_b, path, col_names = T, na = "")
path = paste0(out_path, "4sets_sigGO_withgene_cladespecific_2023-03-19.tsv")

colp = brewer.pal(4, 'Dark2')
# colp2 = c(rep(colp[1], 2), rep(colp[2], 2), colp[3])

N = 3
size_max = ceiling(max(geneCount$n_genes_perclade_perterm)/N)*N
# x = go_anno_gene_sel_cn2 %>% filter(goID =="GO:0051179")

theme_a = theme(axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1),
                axis.text.y = element_text(hjust = 1),
                legend.title = element_text(size = AxisTxFontSizeSize),
                legend.position = "right")

# color of dots in grey ====
# 2.2.1 sort by n_clades, sum_gene ====
y = geneCount_b %>% arrange(n_clades, sum_gene)
y_lab = unique(y$lab)
y$lab = factor(y$lab, y_lab)
y_clade_sp = y

size_max = 3

# 3
p_term3 = y %>% ggplot() +
  geom_point(aes(y = lab, #fct_reorder2(lab, n_genes_perclade_perterm, sum_gene, .desc = F), 
                 x = clade, size = n_genes_perclade_perterm), color = 'grey50') + # color = clade, 
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_x_discrete(name = "# of CNEEs in accelerated state") +
  scale_color_manual(values = colp, guide = "none") +
  scale_size_continuous(name = "Number of genes", range = c(0.2, size_max)) +
  # guides() +
  theme_m + 
  theme_a; p_term3

# use p_term23's order
# x3 = fct_reorder2(geneCount$lab, geneCount$n_genes_perclade_perterm, geneCount$sum_gene, .desc = F)
# x3 = geneCount %>% group_by(ont) %>% mutate(lab = fct_reorder2(lab, n_genes_perclade_perterm, sum_gene, .desc = F))
# x3_2 = x3$lab

pal = brewer.pal(9, 'Reds')
p3 = ggplot( y ) +
  geom_tile(aes(x = 1, 
                y = lab, #fct_relevel(lab, levels(x3_2)), 
                fill = `-log10P`), color = 'white') +
  scale_x_continuous(expand = c(0,0,0.05,0), breaks = 1, labels = expression(-~log[10]~P)) +
  scale_y_discrete(expand = c(0,0,0,0)) +
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_fill_manual(values = pal1) +
  scale_fill_gradient(low = pal[2], high = "#dc2943ff", name = expression(-~log[10]~P)) +
  # labs(x=expression(-~log[10]~P)) +
  theme_m + 
  theme_a; p3

p = cowplot::plot_grid(p3, p_term3, ncol= 2, scale = rep(c(0.9566691, 1), 3), align = 'v'); 
# rel_widths=c(1.3, 1.3, 1, 1.3)) # labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2)

graph_path = paste0(out_path, "perterm_Ngeneperclade_cladespecific_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*1.8*2, height = Width_HalfCol*7.7, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p
dev.off()
openfile(graph_path)

# 2.2.1 sort by n_clades, sum_gene [for terms with genes from all 4 clades] ====
y = geneCount_b %>% filter(n_clades == 4) %>% arrange(n_clades, sum_gene)
y_lab = unique(y$lab)
y$lab = factor(y$lab, y_lab)
y_clade_sp = y

size_max = 3

# 3
p_term3 = y %>% ggplot() +
  geom_point(aes(y = lab, #fct_reorder2(lab, n_genes_perclade_perterm, sum_gene, .desc = F), 
                 x = clade, size = n_genes_perclade_perterm), color = 'grey50') + # color = clade, 
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_x_discrete(name = "# of CNEEs in accelerated state") +
  scale_color_manual(values = colp, guide = "none") +
  scale_size_continuous(name = "Number of genes", range = c(0.2, size_max)) +
  # guides() +
  theme_m + 
  theme_a; p_term3

# use p_term23's order
# x3 = fct_reorder2(geneCount$lab, geneCount$n_genes_perclade_perterm, geneCount$sum_gene, .desc = F)
# x3 = geneCount %>% group_by(ont) %>% mutate(lab = fct_reorder2(lab, n_genes_perclade_perterm, sum_gene, .desc = F))
# x3_2 = x3$lab

pal = brewer.pal(9, 'Reds')
p3 = ggplot( y ) +
  geom_tile(aes(x = 1, 
                y = lab, #fct_relevel(lab, levels(x3_2)), 
                fill = `-log10P`), color = 'white') +
  scale_x_continuous(expand = c(0,0,0.05,0), breaks = 1, labels = expression(-~log[10]~P)) +
  scale_y_discrete(expand = c(0,0,0,0)) +
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_fill_manual(values = pal1) +
  scale_fill_gradient(low = pal[2], high = "#dc2943ff", name = expression(-~log[10]~P)) +
  # labs(x=expression(-~log[10]~P)) +
  theme_m + 
  theme_a; p3

p = cowplot::plot_grid(p3, p_term3, ncol= 2, scale = rep(c(0.9566691, 1), 3), align = 'v'); 
# rel_widths=c(1.3, 1.3, 1, 1.3)) # labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2)

graph_path = paste0(out_path, "perterm_Ngeneperclade_cladespecific_4way_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*3.05, height = Width_HalfCol*1.3, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p
dev.off()
openfile(graph_path)

# 2.2.1 sort by n_clades, sum_gene [for terms with genes from more than 2] ====
y = geneCount_b %>% filter(n_clades >=2) %>% arrange(n_clades, sum_gene)
y_lab = unique(y$lab)
y$lab = factor(y$lab, y_lab)
y_clade_sp = y

size_max = 3

# 3
p_term3 = y %>% ggplot() +
  geom_point(aes(y = lab, #fct_reorder2(lab, n_genes_perclade_perterm, sum_gene, .desc = F), 
                 x = clade, size = n_genes_perclade_perterm), color = 'grey50') + # color = clade, 
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_x_discrete(name = "# of CNEEs in accelerated state") +
  scale_color_manual(values = colp, guide = "none") +
  scale_size_continuous(name = "Number of genes", range = c(0.2, size_max)) +
  # guides() +
  theme_m + 
  theme_a; p_term3

# use p_term23's order
# x3 = fct_reorder2(geneCount$lab, geneCount$n_genes_perclade_perterm, geneCount$sum_gene, .desc = F)
# x3 = geneCount %>% group_by(ont) %>% mutate(lab = fct_reorder2(lab, n_genes_perclade_perterm, sum_gene, .desc = F))
# x3_2 = x3$lab

pal = brewer.pal(9, 'Reds')
p3 = ggplot( y ) +
  geom_tile(aes(x = 1, 
                y = lab, #fct_relevel(lab, levels(x3_2)), 
                fill = `-log10P`), color = 'white') +
  scale_x_continuous(expand = c(0,0,0.05,0), breaks = 1, labels = expression(-~log[10]~P)) +
  scale_y_discrete(expand = c(0,0,0,0)) +
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_fill_manual(values = pal1) +
  scale_fill_gradient(low = pal[2], high = "#dc2943ff", name = expression(-~log[10]~P)) +
  # labs(x=expression(-~log[10]~P)) +
  theme_m + 
  theme_a; p3

p = cowplot::plot_grid(p3, p_term3, ncol= 2, scale = rep(c(0.9566691, 1), 3), align = 'v'); 
# rel_widths=c(1.3, 1.3, 1, 1.3)) # labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2)

graph_path = paste0(out_path, "perterm_Ngeneperclade_cladespecific_gt2way_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*3.05, height = Width_HalfCol*4.5, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p
dev.off()
openfile(graph_path)

# 2.2.1 sort by n_clades, sum_gene [for terms with genes from only 1 clade ] ====
y = geneCount_b %>% filter(n_clades ==1) %>% arrange(n_clades, sum_gene)
y_lab = unique(y$lab)
y$lab = factor(y$lab, y_lab)
y_clade_sp = y

size_max = 3

# 3
p_term3 = y %>% ggplot() +
  geom_point(aes(y = lab, #fct_reorder2(lab, n_genes_perclade_perterm, sum_gene, .desc = F), 
                 x = clade, size = n_genes_perclade_perterm), color = 'grey50') + # color = clade, 
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_x_discrete(name = "# of CNEEs in accelerated state") +
  scale_color_manual(values = colp, guide = "none") +
  scale_size_continuous(name = "Number of genes", range = c(0.2, size_max)) +
  # guides() +
  theme_m + 
  theme_a; p_term3

# use p_term23's order
# x3 = fct_reorder2(geneCount$lab, geneCount$n_genes_perclade_perterm, geneCount$sum_gene, .desc = F)
# x3 = geneCount %>% group_by(ont) %>% mutate(lab = fct_reorder2(lab, n_genes_perclade_perterm, sum_gene, .desc = F))
# x3_2 = x3$lab

pal = brewer.pal(9, 'Reds')
p3 = ggplot( y ) +
  geom_tile(aes(x = 1, 
                y = lab, #fct_relevel(lab, levels(x3_2)), 
                fill = `-log10P`), color = 'white') +
  scale_x_continuous(expand = c(0,0,0.05,0), breaks = 1, labels = expression(-~log[10]~P)) +
  scale_y_discrete(expand = c(0,0,0,0)) +
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_fill_manual(values = pal1) +
  scale_fill_gradient(low = pal[2], high = "#dc2943ff", name = expression(-~log[10]~P)) +
  # labs(x=expression(-~log[10]~P)) +
  theme_m + 
  theme_a; p3

p = cowplot::plot_grid(p3, p_term3, ncol= 2, scale = rep(c(0.9566691, 1), 3), align = 'v'); 
# rel_widths=c(1.3, 1.3, 1, 1.3)) # labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2)

graph_path = paste0(out_path, "perterm_Ngeneperclade_cladespecific_truespecific_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*3.75, height = Width_HalfCol*4.5, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p
dev.off()
openfile(graph_path)



# comparing: clade-specific true specific and 4way ====
y = geneCount_b %>% filter(n_clades == 4) %>% arrange(n_clades, sum_gene)
y_lab = unique(y$lab)
y$lab = factor(y$lab, y_lab)
y_clade_sp4 = y

y = geneCount_b %>% filter(n_clades == 1) %>% arrange(n_clades, sum_gene)
y_lab = unique(y$lab)
y$lab = factor(y$lab, y_lab)
y_clade_sp1 = y

compute_go_sim = sapply(unique(y_clade_sp$ont), function(s){
  tmp_ont = case_when(
    s == "biological_process" ~ "BP",
    s == "molecular_function" ~ "MF",
    s == "cellular_component" ~ "CC" )
  tmp_clade_sp = y_clade_sp4 %>% filter(ont == s)
  tmp_4sets = y_clade_sp1 %>% filter(ont == s)
  
  hsGO <- GOSemSim::godata('org.Hs.eg.db', ont=tmp_ont)
  GOSemSim::mgoSim(tmp_clade_sp$go_id, tmp_4sets$go_id, semData=hsGO, measure="Wang", combine="BMA")
}); compute_go_sim
# biological_process cellular_component molecular_function 
# 0.447              0.556              0.240

# comparing: clade-specific 4way and clade-union ====
y = geneCount_b %>% filter(n_clades == 4) %>% arrange(n_clades, sum_gene)
y_lab = unique(y$lab)
y$lab = factor(y$lab, y_lab)
y_clade_sp4 = y

y_4sets

compute_go_sim = sapply(unique(y_clade_sp$ont), function(s){
  tmp_ont = case_when(
    s == "biological_process" ~ "BP",
    s == "molecular_function" ~ "MF",
    s == "cellular_component" ~ "CC" )
  tmp_clade_sp = y_clade_sp4 %>% filter(ont == s)
  tmp_4sets = y_4sets %>% filter(ont == s)
  
  hsGO <- GOSemSim::godata('org.Hs.eg.db', ont=tmp_ont)
  GOSemSim::mgoSim(tmp_clade_sp$go_id, tmp_4sets$go_id, semData=hsGO, measure="Wang", combine="BMA")
}); compute_go_sim
# biological_process cellular_component molecular_function 
# 0.494              0.628              0.223


# comparing: clade-union vs clade-specific ====
y_clade_sp
y_4sets

compute_go_sim = sapply(unique(y_clade_sp$ont), function(s){
  tmp_ont = case_when(
    s == "biological_process" ~ "BP",
    s == "molecular_function" ~ "MF",
    s == "cellular_component" ~ "CC" )
  tmp_clade_sp = y_clade_sp %>% filter(ont == s)
  tmp_4sets = y_4sets %>% filter(ont == s)
  
  hsGO <- GOSemSim::godata('org.Hs.eg.db', ont=tmp_ont)
  GOSemSim::mgoSim(tmp_clade_sp$go_id, tmp_4sets$go_id, semData=hsGO, measure="Wang", combine="BMA")
})
# target:
# biological_process cellular_component molecular_function 
# 0.595              0.635              0.466 

# control:
# biological_process cellular_component molecular_function 
# 0.639              0.690              0.670 

# comparing: clade-union_4-way vs clade-specific_4way ====
y_4sets_4way = y_4sets %>% filter(n_clades == 4)

y = geneCount_b %>% filter(n_clades == 4) %>% arrange(n_clades, sum_gene)
y_lab = unique(y$lab)
y$lab = factor(y$lab, y_lab)
y_clade_sp4_4way = y

compute_go_sim = sapply(unique(y_clade_sp$ont), function(s){
  tmp_ont = case_when(
    s == "biological_process" ~ "BP",
    s == "molecular_function" ~ "MF",
    s == "cellular_component" ~ "CC" )
  tmp_4sets_con = y_clade_sp4_4way %>% filter(ont == s)
  tmp_4sets = y_4sets_4way %>% filter(ont == s)
  
  hsGO <- GOSemSim::godata('org.Hs.eg.db', ont=tmp_ont)
  GOSemSim::mgoSim(tmp_4sets_con$go_id, tmp_4sets$go_id, semData=hsGO, measure="Wang", combine="BMA")
})
compute_go_sim
# biological_process cellular_component molecular_function 
# 0.539              0.630              0.213 

# comparing: clade-union_4-way vs clade-specific_truesp ====
y_4sets_4way = y_4sets %>% filter(n_clades == 4)

y = geneCount_b %>% filter(n_clades == 1) %>% arrange(n_clades, sum_gene)
y_lab = unique(y$lab)
y$lab = factor(y$lab, y_lab)
y_clade_sp4_1way = y

compute_go_sim = sapply(unique(y_clade_sp$ont), function(s){
  tmp_ont = case_when(
    s == "biological_process" ~ "BP",
    s == "molecular_function" ~ "MF",
    s == "cellular_component" ~ "CC" )
  tmp_4sets_con = y_clade_sp4_1way %>% filter(ont == s)
  tmp_4sets = y_4sets_4way %>% filter(ont == s)
  
  hsGO <- GOSemSim::godata('org.Hs.eg.db', ont=tmp_ont)
  GOSemSim::mgoSim(tmp_4sets_con$go_id, tmp_4sets$go_id, semData=hsGO, measure="Wang", combine="BMA")
})
compute_go_sim
# biological_process cellular_component molecular_function 
#              0.459              0.539              0.419 


# compare target and control: 4sets ====
path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis_2021/05_perm_phyloacc/1_rmPar/postPhyloAcc/go_perms_control/hsa_ensembl_v6v3_2/4sets_sigGO_withgene_control_cladespecific_2023-03-19.tsv")

geneCount_b_con = read_tsv(path)
y = geneCount_b_con %>% arrange(n_clades, sum_gene)
y_lab = unique(y$lab)
y$lab = factor(y$lab, y_lab)
y_4sets_con = y

compute_go_sim = sapply(unique(y_clade_sp$ont), function(s){
  tmp_ont = case_when(
    s == "biological_process" ~ "BP",
    s == "molecular_function" ~ "MF",
    s == "cellular_component" ~ "CC" )
  tmp_4sets_con = y_4sets_con %>% filter(ont == s)
  tmp_4sets = y_4sets %>% filter(ont == s)
  
  hsGO <- GOSemSim::godata('org.Hs.eg.db', ont=tmp_ont)
  GOSemSim::mgoSim(tmp_4sets_con$go_id, tmp_4sets$go_id, semData=hsGO, measure="Wang", combine="BMA")
})
# biological_process cellular_component molecular_function 
# 0.528              0.590              0.479 

# compare 4-way 
y_4sets_con_4way = y_4sets_con %>% filter(n_clades == 4)
y_4sets_4way = y_4sets %>% filter(n_clades == 4)

compute_go_sim = sapply(unique(y_clade_sp$ont), function(s){
  tmp_ont = case_when(
    s == "biological_process" ~ "BP",
    s == "molecular_function" ~ "MF",
    s == "cellular_component" ~ "CC" )
  tmp_4sets_con = y_4sets_con_4way %>% filter(ont == s)
  tmp_4sets = y_4sets_4way %>% filter(ont == s)
  
  hsGO <- GOSemSim::godata('org.Hs.eg.db', ont=tmp_ont)
  GOSemSim::mgoSim(tmp_4sets_con$go_id, tmp_4sets$go_id, semData=hsGO, measure="Wang", combine="BMA")
})
compute_go_sim
# biological_process cellular_component molecular_function 
# 0.445              0.597              0.447

# 4. clade specific - gene lv [with genes] ====
go_anno_splitgene_clade_cladesp = go_anno_splitgene_clade %>% 
  filter( grepl("cladesp_genelv", group)) %>% 
  mutate(lab = fct_drop(lab))
unique(go_anno_splitgene_clade_cladesp$clade)
unique(go_anno_splitgene_clade_cladesp$group)

# v2: num cnee per term per clade
geneCount = go_anno_splitgene_clade_cladesp %>% 
  group_by(go_id, clade) %>% # per term per clade, how many (unique) cnees, (sum of all genes)
  mutate(n_genes_perclade_perterm = length(unique(gene))) %>% 
  group_by(go_id) %>% 
  mutate(n_clades = length(unique(clade)),
         sum_gene = sum(n_genes_perclade_perterm)) 

geneCount_b = geneCount %>% 
  dplyr::select(-gene) %>% distinct() %>% 
  group_by(go_id) %>% 
  mutate(n_clades = length(unique(clade)),
         sum_gene = sum(n_genes_perclade_perterm)) 
dim(geneCount_b) # 4612  17

path = paste0(out_path, "4sets_sigGO_withgene_cladespecific_genelv_", Sys.Date(), ".tsv")
# write_tsv(geneCount_b, path, col_names = T, na = "")
path = paste0(out_path, "4sets_sigGO_withgene_cladespecific_genelv_2023-04-11.tsv")

colp = brewer.pal(4, 'Dark2')
# colp2 = c(rep(colp[1], 2), rep(colp[2], 2), colp[3])

N = 3
size_max = ceiling(max(geneCount$n_genes_perclade_perterm)/N)*N
# x = go_anno_gene_sel_cn2 %>% filter(goID =="GO:0051179")

theme_a = theme(axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1),
                axis.text.y = element_text(hjust = 1),
                legend.title = element_text(size = AxisTxFontSizeSize),
                legend.position = "right")

# color of dots in grey ====
# 2.2.1 sort by n_clades, sum_gene ====
y = geneCount_b %>% arrange(n_clades, sum_gene)
y_lab = unique(y$lab)
y$lab = factor(y$lab, y_lab)
y_clade_sp_genelv = y

size_max = 3

# 3
p_term3 = y %>% ggplot() +
  geom_point(aes(y = lab, #fct_reorder2(lab, n_genes_perclade_perterm, sum_gene, .desc = F), 
                 x = clade, size = n_genes_perclade_perterm), color = 'grey50') + # color = clade, 
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_x_discrete(name = "# of CNEEs in accelerated state") +
  scale_color_manual(values = colp, guide = "none") +
  scale_size_continuous(name = "Number of genes", range = c(0.2, size_max)) +
  # guides() +
  theme_m + 
  theme_a; p_term3

# use p_term23's order
# x3 = fct_reorder2(geneCount$lab, geneCount$n_genes_perclade_perterm, geneCount$sum_gene, .desc = F)
# x3 = geneCount %>% group_by(ont) %>% mutate(lab = fct_reorder2(lab, n_genes_perclade_perterm, sum_gene, .desc = F))
# x3_2 = x3$lab

pal = brewer.pal(9, 'Reds')
p3 = ggplot( y ) +
  geom_tile(aes(x = 1, 
                y = lab, #fct_relevel(lab, levels(x3_2)), 
                fill = `-log10P`), color = 'white') +
  scale_x_continuous(expand = c(0,0,0.05,0), breaks = 1, labels = expression(-~log[10]~P)) +
  scale_y_discrete(expand = c(0,0,0,0)) +
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_fill_manual(values = pal1) +
  scale_fill_gradient(low = pal[2], high = "#dc2943ff", name = expression(-~log[10]~P)) +
  # labs(x=expression(-~log[10]~P)) +
  theme_m + 
  theme_a; p3

p = cowplot::plot_grid(p3, p_term3, ncol= 2, scale = rep(c(0.9566691, 1), 3), align = 'v'); 
# rel_widths=c(1.3, 1.3, 1, 1.3)) # labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2)

graph_path = paste0(out_path, "perterm_Ngeneperclade_cladespecific_genelv_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*1.8*2, height = Width_HalfCol*7.7, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p
dev.off()
openfile(graph_path)

# 2.2.1 sort by n_clades, sum_gene [for terms with genes from more than 2] ====
y = geneCount_b %>% filter(n_clades >=2) %>% arrange(n_clades, sum_gene)
y_lab = unique(y$lab)
y$lab = factor(y$lab, y_lab)
y_clade_sp = y

size_max = 3

# 3
p_term3 = y %>% ggplot() +
  geom_point(aes(y = lab, #fct_reorder2(lab, n_genes_perclade_perterm, sum_gene, .desc = F), 
                 x = clade, size = n_genes_perclade_perterm), color = 'grey50') + # color = clade, 
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_x_discrete(name = "# of CNEEs in accelerated state") +
  scale_color_manual(values = colp, guide = "none") +
  scale_size_continuous(name = "Number of genes", range = c(0.2, size_max)) +
  # guides() +
  theme_m + 
  theme_a; p_term3

# use p_term23's order
# x3 = fct_reorder2(geneCount$lab, geneCount$n_genes_perclade_perterm, geneCount$sum_gene, .desc = F)
# x3 = geneCount %>% group_by(ont) %>% mutate(lab = fct_reorder2(lab, n_genes_perclade_perterm, sum_gene, .desc = F))
# x3_2 = x3$lab

pal = brewer.pal(9, 'Reds')
p3 = ggplot( y ) +
  geom_tile(aes(x = 1, 
                y = lab, #fct_relevel(lab, levels(x3_2)), 
                fill = `-log10P`), color = 'white') +
  scale_x_continuous(expand = c(0,0,0.05,0), breaks = 1, labels = expression(-~log[10]~P)) +
  scale_y_discrete(expand = c(0,0,0,0)) +
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_fill_manual(values = pal1) +
  scale_fill_gradient(low = pal[2], high = "#dc2943ff", name = expression(-~log[10]~P)) +
  # labs(x=expression(-~log[10]~P)) +
  theme_m + 
  theme_a; p3

p = cowplot::plot_grid(p3, p_term3, ncol= 2, scale = rep(c(0.9566691, 1), 3), align = 'v'); 
# rel_widths=c(1.3, 1.3, 1, 1.3)) # labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2)

graph_path = paste0(out_path, "perterm_Ngeneperclade_cladespecific_genelv_gt2way_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*5, height = Width_HalfCol*20, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p
dev.off()
openfile(graph_path)


# 2.2.1 sort by n_clades, sum_gene [for terms with genes from only 1 clade ] ====
y = geneCount_b %>% filter(n_clades ==1) %>% arrange(n_clades, sum_gene)
y_lab = unique(y$lab)
y$lab = factor(y$lab, y_lab)
y_clade_sp = y

size_max = 3

# 3
p_term3 = y %>% ggplot() +
  geom_point(aes(y = lab, #fct_reorder2(lab, n_genes_perclade_perterm, sum_gene, .desc = F), 
                 x = clade, size = n_genes_perclade_perterm), color = 'grey50') + # color = clade, 
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_x_discrete(name = "# of CNEEs in accelerated state") +
  scale_color_manual(values = colp, guide = "none") +
  scale_size_continuous(name = "Number of genes", range = c(0.2, size_max)) +
  # guides() +
  theme_m + 
  theme_a; p_term3

# use p_term23's order
# x3 = fct_reorder2(geneCount$lab, geneCount$n_genes_perclade_perterm, geneCount$sum_gene, .desc = F)
# x3 = geneCount %>% group_by(ont) %>% mutate(lab = fct_reorder2(lab, n_genes_perclade_perterm, sum_gene, .desc = F))
# x3_2 = x3$lab

pal = brewer.pal(9, 'Reds')
p3 = ggplot( y ) +
  geom_tile(aes(x = 1, 
                y = lab, #fct_relevel(lab, levels(x3_2)), 
                fill = `-log10P`), color = 'white') +
  scale_x_continuous(expand = c(0,0,0.05,0), breaks = 1, labels = expression(-~log[10]~P)) +
  scale_y_discrete(expand = c(0,0,0,0)) +
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_fill_manual(values = pal1) +
  scale_fill_gradient(low = pal[2], high = "#dc2943ff", name = expression(-~log[10]~P)) +
  # labs(x=expression(-~log[10]~P)) +
  theme_m + 
  theme_a; p3

p = cowplot::plot_grid(p3, p_term3, ncol= 2, scale = rep(c(0.9566691, 1), 3), align = 'v'); 
# rel_widths=c(1.3, 1.3, 1, 1.3)) # labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2)

graph_path = paste0(out_path, "perterm_Ngeneperclade_cladespecific_genelv_truespecific_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*5, height = Width_HalfCol*30, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p
dev.off()
openfile(graph_path)

## C.  get gene symbol ====
## C1. get ori gene symbol 
# get entrez to gene symbol (from the ori ncbi galgal6 prepared for permutation)
target_files = list.files(project_path, ".counts.bed$", full.names = T); target_files
get_galgal_func = function(s){
  # s = target_files[1]
  print(paste0("current: ", s))
  tmp_name = basename(s)
  tmp_name = gsub("cnee_goperms_(.*)\\.counts\\.bed", "\\1", tmp_name )
  print(paste0("current: ", tmp_name))
  
  # load in background data and clean up NCBI IDs
  bg <- read_tsv(s, col_names = F, show_col_types = FALSE) %>%
    dplyr::rename(chr = X1, start = X2, end = X3, gene = X4, combo = X5, biotype = X6, accel = X7, total = X8) %>%
    separate(combo, into = c(NA, "ncbi"), sep = "GeneID:", remove = T) %>%
    dplyr::select(gene, ncbi)
}
ori_geneIDs0 = lapply(target_files, get_galgal_func)
ori_geneIDs = ori_geneIDs0 %>% bind_rows()

## C2. split geneID column
# 0926 modify:
# no need to use 'bitr' to get symbol, we can simply use the ori file to get the ori galgal6 gene symbol

go_anno_splitgene = go_anno %>% 
  mutate(geneIDsplit = strsplit(geneID, split = "/")) %>% 
  unnest(geneIDsplit) %>% 
  dplyr::select(-geneID) %>% dplyr::rename(gene = geneIDsplit)

go_anno_gene = go_anno_splitgene %>% 
  left_join(ori_geneIDs, by = c("gene")) %>% 
  distinct() %>% 
  dplyr::select(group:bg_total, gene, ncbi, goTerm_def)
dim(go_anno_gene) #  3312   13
sum(is.na(go_anno_gene$gene)) # 0

path = paste0(out_path, "hsa_go_anno_obs_genes_", Sys.Date(), ".tsv")
# write_tsv(go_anno_gene, path, col_names = T, na = "")
path = paste0(out_path, "hsa_go_anno_obs_genes_2022-11-17.tsv")
# go_anno_gene = read_tsv(path)
# hsa = go_anno_gene


# D. load cnee and genes with acc. cnees around (100 kb) ====
target_file = list.files(paste0(project_path, "../linear/"), 
                         "cnee_ncbigene100kb_acc_intersect_FullTreeExp_\\w+.bed", 
                         full.names = T); target_file
acccnees_100kbgg_withcnees = lapply(target_file, function(s){
  # s = target_file[1]
  tmp_name = basename(s)
  tmp_name = gsub("cnee_ncbigene100kb_acc_intersect_FullTreeExp_", "", tmp_name)
  tmp_name = gsub(".bed", "", tmp_name)
  print(paste0("current: ", tmp_name))
  x = read_tsv(s, col_names = F) %>% 
    filter(X6 == "protein_coding") %>%
    dplyr::select(X4, X5, X12, X7, X8) %>%
    dplyr::rename("gene" = "X4",  "combo" = "X5", "id" = "X12", 
                  "accel" = "X7", "total" = "X8") %>%
    separate(combo, into = c(NA, "ncbi"), sep = "GeneID:", remove = T) %>%
    bind_cols('group' = tmp_name) %>% 
    dplyr::select(group, id, gene, accel, total) # ncbi isn't better
}) %>% bind_rows()
dim(acccnees_100kbgg_withcnees) # 18936     5
unique(acccnees_100kbgg_withcnees$group)
# [1] "2wayconvergent"        "3wayconvergent"        "4sets"                 "genelv_conv_2way"      "genelv_conv_3way"     
# [6] "genelv_conv_4way"      "honeyeaters"           "hummingbirds"          "parrots"               "sunbirds"
head(acccnees_100kbgg_withcnees)
sum(is.na(acccnees_100kbgg_withcnees$gene)) # 0
# sum(is.na(acccnees_100kbgg_withcnees$ncbi)) # 0
sum(is.na(acccnees_100kbgg_withcnees$id)) # 0

go_obsGene_cnee = go_anno_gene %>% 
  left_join(acccnees_100kbgg_withcnees, by = c("gene", "group")) # better
  # left_join(acccnees_100kbgg_withcnees, by = c("geneID" = "ncbi", "group"))
dim(go_anno_gene) # 3312   13
dim(go_obsGene_cnee) # 5525   16
head(go_obsGene_cnee)
# target_in: number of genes observed in GOterm
# accel: number of acc. cnees for that gene in this group
# total: total of acc. cnees for that gene
# gene is the original galgal annotation from the GO count file

# check
sum(is.na(go_obsGene_cnee$gene)) # 0
# sum(is.na(go_obsGene_cnee$gene.y)) # 2
sum(is.na(go_obsGene_cnee$id)) # left_join using gene: 3


x = go_obsGene_cnee %>% filter(is.na(gene)) # 0
y = go_obsGene_cnee %>% filter(is.na(id)) # GRXCR2
# unique(y$SYMBOL) # GRXCR2
unique(y$gene) # GRXCR2, not sure why

path = paste0(out_path, "Final_GO_combined_sig_obsgene_cnees_", Sys.Date(), ".tsv")
# write_tsv(go_obsGene_cnee, path, col_names = T, na = "")
# path = paste0(out_path, "Final_GO_combined_sig_obsgene_cnees_2022-11-17.tsv")
# go_obsGene_cnee = read_tsv(path)

##### investigating =====

acccnees_100kbgg_withcnees_4 = acccnees_100kbgg_withcnees %>% 
  filter(group %in% c("honeyeaters", "hummingbirds", "parrots", "sunbirds")) %>% 
  dplyr::select(gene, id, group) %>% distinct()

# check
x = go_obsGene_cnee %>% filter(group == "2wayconvergent", go_id == "GO:0001934")
y = acccnees_100kbgg_withcnees %>% filter(group == "2wayconvergent", gene %in% x$gene)
# ok, x's cnees are identical as y's cnees

###### genelv_conv_4way ====
unique(go_obsGene_cnee$group)
tmp_set = 'genelv_conv_4way'
go_anno_gene_sel <- go_anno_gene %>% filter(group == tmp_set) 

go_anno_gene_sel_cn = go_anno_gene_sel %>% 
  # this left_join is necessary, to get each clade's cnees
  left_join(acccnees_100kbgg_withcnees_4, by = c("gene")) %>%
  # filter( !is.na(id)) %>% # not necessary because no na
  group_by(go_id) %>%  # per go_id
  mutate(n_clades_per_term = length(unique(group.y))) %>% 
  group_by(go_id, gene) %>%  # per gene
  mutate(n_clades_per_gene = length(unique(group.y))) %>% 
  group_by(go_id, gene, id) %>% # per cnee
  mutate(n_clades_per_cnee = length(unique(group.y))) %>% 
  mutate(lab = paste0(go_id, ": ", goTerm)) # %>%
  # mutate(lab = fct_reorder(lab,  total, .desc = T)) 

length(unique(go_anno_gene_sel_cn$gene)) # 4
length(unique(go_anno_gene_sel_cn$id)) # 29

# go_anno_gene_sel_cn$lab = factor(go_anno_gene_sel_cn$lab, 
#                                   levels = c("GO:0005488\nbinding", 
# "GO:0043226\norganelle", "GO:0043229\nintracellular organelle",
# "GO:0043228\nnon-membrane-bounded organelle", 
# "GO:0043232\nintracellular non-membrane-bounded organelle"))

# go_anno_gene_sel_cn = go_anno_gene_sel_cn %>%
# mutate(lab = fct_reorder(lab,  n_clades_per_term, .desc = T)) # doesnt work here, dont know why
# go_anno_gene_sel_cn$lab = factor(go_anno_gene_sel_cn$lab, 
#                                   levels = sort(unique(go_anno_gene_sel_cn$lab))) # same as defualt

go_anno_gene_sel_cn$group.y = factor(go_anno_gene_sel_cn$group.y, 
                                      levels = c("hummingbirds", "parrots", 
                                                 "honeyeaters", "sunbirds"))

# v2: num cnee per term per clade
go_anno_gene_sel_cn2 = go_anno_gene_sel_cn %>% 
  group_by(go_id, group.y) %>% # per term per clade, how many (unique) cnees, (sum of all genes)
  mutate(n_cnees_perclade_perterm = length(unique(id)))

length(unique(go_anno_gene_sel_cn2$ont)) # 2

colp = brewer.pal(4, 'Dark2')
# colp2 = c(rep(colp[1], 2), rep(colp[2], 2), colp[3])

N = 4
size_max = ceiling(max(go_anno_gene_sel_cn2$n_cnees_perclade_perterm)/N)*N
x = go_anno_gene_sel_cn2 %>% filter(go_id =="GO:0051179")

p_term2 = go_anno_gene_sel_cn2 %>% ggplot() +
  geom_point(aes(y = lab, x = group.y, color = group.y, size = n_cnees_perclade_perterm)) +
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_x_discrete(name = "# of CNEEs in accelerated state") +
  scale_color_manual(values = colp, guide = "none") +
  scale_size_continuous(name = "Number of accelerated CNEEs", range = c(0.2, 4), limits = c(1, size_max)) +
  # guides() +
  theme_m + 
  theme(#plot.margin = unit(rep(4,4), "pt"), # trbl
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1),
    axis.text.y = element_text(hjust = 1),
    legend.position = "bottom"); p_term2

graph_path = paste0(out_path, "GO_diff_convergent_lv_", tmp_set, "_perterm_Ncneeperclade_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*1.7, height = Width_HalfCol*0.75, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p_term2
dev.off()
openfile(graph_path)

# v3: num cnee per term per clade, gene lv
go_anno_gene_sel_cn2 = go_anno_gene_sel_cn %>% 
  group_by(go_id, group.y) %>% # per term per clade, how many (unique) cnees, (sum of all genes)
  mutate(n_cnees_perclade_perterm = length(unique(id))) %>% 
  group_by(go_id, group.y, gene) %>% # per term per clade, how many (unique) cnees, per gene
  mutate(n_cnees_perclade_perterm_pergene = length(unique(id)))

colp = brewer.pal(4, 'Dark2')
# colp2 = c(rep(colp[1], 2), rep(colp[2], 2), colp[3])

N = 4
size_max = ceiling(max(go_anno_gene_sel_cn2$n_cnees_perclade_perterm_pergene)/N)*N
x = go_anno_gene_sel_cn2 %>% filter(go_id =="GO:0051179")

p_term3 = go_anno_gene_sel_cn2 %>% 
  # filter(ont == "BP") %>% 
  # filter(n_clades_per_gene >= 2) %>% 
  ggplot() +
  geom_point(aes(y = gene, x = group.y, color = group.y, size = n_cnees_perclade_perterm_pergene)) +
  facet_grid(ont+lab~.  , scales = 'free', space = 'free') +
  # facet_wrap(lab~.  , scales = 'free', ncol = 3) +
  # scale_x_discrete(name = "# of CNEEs in accelerated state") +
  scale_y_discrete(limits = rev) +
  scale_color_manual(values = colp) + # , guide = "none"
  scale_size_continuous(name = "Number of accelerated CNEEs", range = c(0.2, 4), limits = c(1, size_max)) +
  guides(color = guide_legend(nrow = 4)) +
  theme_m + 
  theme(#plot.margin = unit(rep(4,4), "pt"), # trbl
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1),
    axis.text.y = element_text(hjust = 1),
    strip.text.y = element_text(angle = 0, hjust = 0),
    legend.position = "bottom"); p_term3

graph_path = paste0(out_path, "GO_diff_convergent_lv_", tmp_set, "_perterm_Ncneeperclade_genelv_bp_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*3, height = Width_HalfCol*1, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p_term3
dev.off()
openfile(graph_path)

###### genelv_conv_3way ====
unique(go_obsGene_cnee$group)
tmp_set = 'genelv_conv_3way'
go_gene_list_sel2 <- go_obsGene_cnee %>% filter(group == tmp_set) 
sum(is.na(go_gene_list_sel2$id)) # 0

go_anno_gene_sel_cn = go_gene_list_sel2 %>% 
  left_join(acccnees_100kbgg_withcnees_4, by = c( "gene", "id")) %>% # this is to add origin's
  # filter( !is.na(id)) %>% 
  group_by(go_id) %>%  # per go_id
  mutate(n_clades_per_term = length(unique(group.y))) %>% 
  group_by(go_id, gene) %>%  # per gene
  mutate(n_clades_per_gene = length(unique(group.y))) %>% 
  group_by(go_id, gene, id) %>% # per cnee
  mutate(n_clades_per_cnee = length(unique(group.y))) %>% 
  mutate(lab = paste0(go_id, ": ", goTerm))  %>%
  mutate(lab = fct_reorder(lab,  total, .desc = T)) 

length(unique(go_anno_gene_sel_cn$gene)) # 57
length(unique(go_anno_gene_sel_cn$id)) # 169


# v2: num cnee per term per clade
go_anno_gene_sel_cn2 = go_anno_gene_sel_cn %>% 
  group_by(go_id, group.y) %>% # per term per clade, how many (unique) cnees, (sum of all genes)
  mutate(n_cnees_perclade_perterm = length(unique(id)))

colp = brewer.pal(4, 'Dark2')
# colp2 = c(rep(colp[1], 2), rep(colp[2], 2), colp[3])

N = 4
size_max = ceiling(max(go_anno_gene_sel_cn2$n_cnees_perclade_perterm)/N)*N
x = go_anno_gene_sel_cn2 %>% filter(go_id =="GO:0051179")

p_term2 = go_anno_gene_sel_cn2 %>% ggplot() +
  geom_point(aes(y = lab, x = group.y, color = group.y, size = n_cnees_perclade_perterm)) +
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_x_discrete(name = "# of CNEEs in accelerated state") +
  scale_color_manual(values = colp, guide = "none") +
  scale_size_continuous(name = "Number of accelerated CNEEs", range = c(0.2, 4), limits = c(1, size_max)) +
  # guides() +
  theme_m + 
  theme(#plot.margin = unit(rep(4,4), "pt"), # trbl
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1),
    axis.text.y = element_text(hjust = 1),
    legend.position = "bottom"); p_term2

graph_path = paste0(out_path, "GO_diff_convergent_lv_", tmp_set, "_perterm_ncneeperclade_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*1.75, height = Width_HalfCol*1.5, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p_term2
dev.off()
openfile(graph_path)

# v3: num cnee per term per clade, gene lv
go_anno_gene_sel_cn2 = go_anno_gene_sel_cn %>% 
  group_by(go_id, group.y) %>% # per term per clade, how many (unique) cnees, (sum of all genes)
  mutate(n_cnees_perclade_perterm = length(unique(id))) %>% 
  group_by(go_id, group.y, gene) %>% # per term per clade, how many (unique) cnees, per gene
  mutate(n_cnees_perclade_perterm_pergene = length(unique(id)))

colp = brewer.pal(4, 'Dark2')
# colp2 = c(rep(colp[1], 2), rep(colp[2], 2), colp[3])

N = 4
size_max = ceiling(max(go_anno_gene_sel_cn2$n_cnees_perclade_perterm_pergene)/N)*N
x = go_anno_gene_sel_cn2 %>% filter(go_id =="GO:0051179")

p_term3 = go_anno_gene_sel_cn2 %>% 
  # filter(ont == "BP") %>% 
  # filter(n_clades_per_gene >= 2) %>% 
  ggplot() +
  geom_point(aes(y = gene, x = group.y, color = group.y, size = n_cnees_perclade_perterm_pergene)) +
  facet_grid(ont+lab~.  , scales = 'free', space = 'free') +
  # facet_wrap(lab~.  , scales = 'free', ncol = 3) +
  # scale_x_discrete(name = "# of CNEEs in accelerated state") +
  scale_y_discrete(limits = rev) +
  scale_color_manual(values = colp) + # , guide = "none"
  scale_size_continuous(name = "Number of accelerated CNEEs", range = c(0.2, 4), limits = c(1, size_max)) +
  guides(color = guide_legend(nrow = 4)) +
  theme_m + 
  theme(#plot.margin = unit(rep(4,4), "pt"), # trbl
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1),
    axis.text.y = element_text(hjust = 1),
    strip.text.y = element_text(angle = 0, hjust = 0),
    legend.position = "bottom"); p_term3

graph_path = paste0(out_path, "GO_diff_convergent_lv_", tmp_set, "_perterm_ncneeperclade_genelv_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2.5, height = Width_HalfCol*4.5, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p_term3
dev.off()
openfile(graph_path)


###### genelv_conv_2way ====
unique(go_obsGene_cnee$group)
tmp_set = 'genelv_conv_2way'
go_gene_list_sel2 <- go_obsGene_cnee %>% filter(group == tmp_set) 
sum(is.na(go_gene_list_sel2$id)) # 0

go_anno_gene_sel_cn = go_gene_list_sel2 %>% 
  left_join(acccnees_100kbgg_withcnees_4, by = c( "gene", "id")) %>% # this is to add origin's
  # filter( !is.na(id)) %>% 
  group_by(go_id) %>%  # per go_id
  mutate(n_clades_per_term = length(unique(group.y))) %>% 
  group_by(go_id, gene) %>%  # per gene
  mutate(n_clades_per_gene = length(unique(group.y))) %>% 
  group_by(go_id, gene, id) %>% # per cnee
  mutate(n_clades_per_cnee = length(unique(group.y))) %>% 
  mutate(lab = paste0(go_id, ": ", goTerm))  %>%
  mutate(lab = fct_reorder(lab,  total, .desc = T)) 

length(unique(go_anno_gene_sel_cn$gene)) # 223
length(unique(go_anno_gene_sel_cn$id)) # 341

# go_anno_gene_sel_cn$lab = factor(go_anno_gene_sel_cn$lab, 
#                                   levels = c("GO:0005488\nbinding", "GO:0043226\norganelle", "GO:0043229\nintracellular organelle",
#                                              "GO:0043228\nnon-membrane-bounded organelle", 
#                                              "GO:0043232\nintracellular non-membrane-bounded organelle"))

# go_anno_gene_sel_cn = go_anno_gene_sel_cn %>%
# mutate(lab = fct_reorder(lab,  n_clades_per_term, .desc = T)) # doesnt work here, dont know why
# go_anno_gene_sel_cn$lab = factor(go_anno_gene_sel_cn$lab, 
#                                   levels = sort(unique(go_anno_gene_sel_cn$lab))) # same as defualt

go_anno_gene_sel_cn$group.y = factor(go_anno_gene_sel_cn$group.y, 
                                      levels = c("hummingbirds", "parrots", "honeyeaters", "sunbirds"))


# v2: num cnee per term per clade
go_anno_gene_sel_cn2 = go_anno_gene_sel_cn %>% 
  group_by(go_id, group.y) %>% # per term per clade, how many (unique) cnees, (sum of all genes)
  mutate(n_cnees_perclade_perterm = length(unique(id)))

colp = brewer.pal(4, 'Dark2')
# colp2 = c(rep(colp[1], 2), rep(colp[2], 2), colp[3])

N = 4
size_max = ceiling(max(go_anno_gene_sel_cn2$n_cnees_perclade_perterm)/N)*N
x = go_anno_gene_sel_cn2 %>% filter(go_id =="GO:0051179")

p_term2 = go_anno_gene_sel_cn2 %>% ggplot() +
  geom_point(aes(y = lab, x = group.y, color = group.y, size = n_cnees_perclade_perterm)) +
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_x_discrete(name = "# of CNEEs in accelerated state") +
  scale_color_manual(values = colp, guide = "none") +
  scale_size_continuous(name = "Number of accelerated CNEEs", range = c(0.2, 4), limits = c(1, size_max)) +
  # guides() +
  theme_m + 
  theme(#plot.margin = unit(rep(4,4), "pt"), # trbl
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1),
    axis.text.y = element_text(hjust = 1),
    legend.position = "bottom"); p_term2

graph_path = paste0(out_path, "GO_diff_convergent_lv_", tmp_set, "_perterm_ncneeperclade_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*1.75, height = Width_HalfCol*2.5, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p_term2
dev.off()
openfile(graph_path)

# v3: num cnee per term per clade, gene lv
go_anno_gene_sel_cn2 = go_anno_gene_sel_cn %>% 
  group_by(go_id, group.y) %>% # per term per clade, how many (unique) cnees, (sum of all genes)
  mutate(n_cnees_perclade_perterm = length(unique(id))) %>% 
  group_by(go_id, group.y, gene) %>% # per term per clade, how many (unique) cnees, per gene
  mutate(n_cnees_perclade_perterm_pergene = length(unique(id)))

colp = brewer.pal(4, 'Dark2')
# colp2 = c(rep(colp[1], 2), rep(colp[2], 2), colp[3])

N = 4
size_max = ceiling(max(go_anno_gene_sel_cn2$n_cnees_perclade_perterm_pergene)/N)*N
x = go_anno_gene_sel_cn2 %>% filter(go_id =="GO:0051179")

p_term3 = go_anno_gene_sel_cn2 %>% 
  # filter(ont == "BP") %>% 
  # filter(n_clades_per_gene >= 2) %>% 
  ggplot() +
  geom_point(aes(y = gene, x = group.y, color = group.y, size = n_cnees_perclade_perterm_pergene)) +
  facet_grid(ont+lab~.  , scales = 'free', space = 'free') +
  # facet_wrap(lab~.  , scales = 'free', ncol = 3) +
  # scale_x_discrete(name = "# of CNEEs in accelerated state") +
  scale_y_discrete(limits = rev) +
  scale_color_manual(values = colp) + # , guide = "none"
  scale_size_continuous(name = "Number of accelerated CNEEs", range = c(0.2, 4), limits = c(1, size_max)) +
  guides(color = guide_legend(nrow = 4)) +
  theme_m + 
  theme(#plot.margin = unit(rep(4,4), "pt"), # trbl
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1),
    axis.text.y = element_text(hjust = 1),
    strip.text.y = element_text(angle = 0, hjust = 0),
    legend.position = "bottom"); p_term3

graph_path = paste0(out_path, "GO_diff_convergent_lv_", tmp_set, "_perterm_ncneeperclade_genelv_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*3.5, height = Width_HalfCol*11, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p_term3
dev.off()
openfile(graph_path)

###### 3wayconvergent ====
unique(go_obsGene_cnee$group)
tmp_set = '3wayconvergent'
go_gene_list_sel2 <- go_obsGene_cnee %>% filter(group == tmp_set) 
sum(is.na(go_gene_list_sel2$id)) # 0

go_anno_gene_sel_cn = go_gene_list_sel2 %>% 
  left_join(acccnees_100kbgg_withcnees_4, by = c( "gene", "id")) %>% # this is to add origin's
  # filter( !is.na(id)) %>% 
  group_by(go_id) %>%  # per go_id
  mutate(n_clades_per_term = length(unique(group.y))) %>% 
  group_by(go_id, gene) %>%  # per gene
  mutate(n_clades_per_gene = length(unique(group.y))) %>% 
  group_by(go_id, gene, id) %>% # per cnee
  mutate(n_clades_per_cnee = length(unique(group.y))) %>% 
  mutate(lab = paste0(go_id, ": ", goTerm))  %>%
  mutate(lab = fct_reorder(lab,  total, .desc = T)) 

length(unique(go_anno_gene_sel_cn$gene)) # 15
length(unique(go_anno_gene_sel_cn$id)) # 10

# go_anno_gene_sel_cn$lab = factor(go_anno_gene_sel_cn$lab, 
#                                   levels = c("GO:0005488\nbinding", "GO:0043226\norganelle", "GO:0043229\nintracellular organelle",
#                                              "GO:0043228\nnon-membrane-bounded organelle", 
#                                              "GO:0043232\nintracellular non-membrane-bounded organelle"))

# go_anno_gene_sel_cn = go_anno_gene_sel_cn %>%
# mutate(lab = fct_reorder(lab,  n_clades_per_term, .desc = T)) # doesnt work here, dont know why
# go_anno_gene_sel_cn$lab = factor(go_anno_gene_sel_cn$lab, 
#                                   levels = sort(unique(go_anno_gene_sel_cn$lab))) # same as defualt

go_anno_gene_sel_cn$group.y = factor(go_anno_gene_sel_cn$group.y, 
                                     levels = c("hummingbirds", "parrots", "honeyeaters", "sunbirds"))

# v2: num cnee per term per clade
go_anno_gene_sel_cn2 = go_anno_gene_sel_cn %>% 
  group_by(go_id, group.y) %>% # per term per clade, how many (unique) cnees, (sum of all genes)
  mutate(n_cnees_perclade_perterm = length(unique(id)))

colp = brewer.pal(4, 'Dark2')
# colp2 = c(rep(colp[1], 2), rep(colp[2], 2), colp[3])

N = 4
size_max = ceiling(max(go_anno_gene_sel_cn2$n_cnees_perclade_perterm)/N)*N
x = go_anno_gene_sel_cn2 %>% filter(go_id =="GO:0051179")

p_term2 = go_anno_gene_sel_cn2 %>% ggplot() +
  geom_point(aes(y = lab, x = group.y, color = group.y, size = n_cnees_perclade_perterm)) +
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_x_discrete(name = "# of CNEEs in accelerated state") +
  scale_color_manual(values = colp, guide = "none") +
  scale_size_continuous(name = "Number of accelerated CNEEs", range = c(0.2, 4), limits = c(1, size_max)) +
  # guides() +
  theme_m + 
  theme(#plot.margin = unit(rep(4,4), "pt"), # trbl
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1),
    axis.text.y = element_text(hjust = 1),
    legend.position = "bottom"); p_term2

graph_path = paste0(out_path, "GO_diff_convergent_lv_", tmp_set, "_perterm_ncneeperclade_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*2.25, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p_term2
dev.off()
openfile(graph_path)

# v3: num cnee per term per clade, gene lv
go_anno_gene_sel_cn2 = go_anno_gene_sel_cn %>% 
  group_by(go_id, group.y) %>% # per term per clade, how many (unique) cnees, (sum of all genes)
  mutate(n_cnees_perclade_perterm = length(unique(id))) %>% 
  group_by(go_id, group.y, gene) %>% # per term per clade, how many (unique) cnees, per gene
  mutate(n_cnees_perclade_perterm_pergene = length(unique(id)))

colp = brewer.pal(4, 'Dark2')
# colp2 = c(rep(colp[1], 2), rep(colp[2], 2), colp[3])

N = 4
size_max = ceiling(max(go_anno_gene_sel_cn2$n_cnees_perclade_perterm_pergene)/N)*N
x = go_anno_gene_sel_cn2 %>% filter(go_id =="GO:0051179")

p_term3 = go_anno_gene_sel_cn2 %>% 
  # filter(ont == "BP") %>% 
  # filter(n_clades_per_gene >= 2) %>% 
  ggplot() +
  geom_point(aes(y = lab, x = group.y, color = group.y, size = n_cnees_perclade_perterm_pergene)) +
  facet_grid(ont+gene~.  , scales = 'free', space = 'free') +
  # facet_wrap(lab~.  , scales = 'free', ncol = 3) +
  # scale_x_discrete(name = "# of CNEEs in accelerated state") +
  scale_y_discrete(limits = rev) +
  scale_color_manual(values = colp) + # , guide = "none"
  scale_size_continuous(name = "Number of accelerated CNEEs", range = c(0.2, 4), limits = c(1, size_max)) +
  guides(color = guide_legend(nrow = 4)) +
  theme_m + 
  theme(#plot.margin = unit(rep(4,4), "pt"), # trbl
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1),
    axis.text.y = element_text(hjust = 1),
    strip.text.y = element_text(angle = 0, hjust = 0),
    legend.position = "bottom"); p_term3

graph_path = paste0(out_path, "GO_diff_convergent_lv_", tmp_set, "_perterm_ncneeperclade_genelv_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2.5, height = Width_HalfCol*2.25, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p_term3
dev.off()
openfile(graph_path)

###### 2wayconvergent ====
unique(go_obsGene_cnee$group)
tmp_set = '2wayconvergent'
go_gene_list_sel2 <- go_obsGene_cnee %>% filter(group == tmp_set) 
sum(is.na(go_gene_list_sel2$id)) # 0

go_anno_gene_sel_cn = go_gene_list_sel2 %>% 
  left_join(acccnees_100kbgg_withcnees_4, by = c( "gene", "id")) %>% # this is to add origin's
  # filter( !is.na(id)) %>% 
  group_by(go_id) %>%  # per go_id
  mutate(n_clades_per_term = length(unique(group.y))) %>% 
  group_by(go_id, gene) %>%  # per gene
  mutate(n_clades_per_gene = length(unique(group.y))) %>% 
  group_by(go_id, gene, id) %>% # per cnee
  mutate(n_clades_per_cnee = length(unique(group.y))) %>% 
  mutate(lab = paste0(go_id, ": ", goTerm))  %>%
  mutate(lab = fct_reorder(lab,  total, .desc = T)) 

length(unique(go_anno_gene_sel_cn$gene)) # 70
length(unique(go_anno_gene_sel_cn$id)) # 67

# go_anno_gene_sel_cn$lab = factor(go_anno_gene_sel_cn$lab, 
#                                   levels = c("GO:0005488\nbinding", "GO:0043226\norganelle", "GO:0043229\nintracellular organelle",
#                                              "GO:0043228\nnon-membrane-bounded organelle", 
#                                              "GO:0043232\nintracellular non-membrane-bounded organelle"))

# go_anno_gene_sel_cn = go_anno_gene_sel_cn %>%
# mutate(lab = fct_reorder(lab,  n_clades_per_term, .desc = T)) # doesnt work here, dont know why
# go_anno_gene_sel_cn$lab = factor(go_anno_gene_sel_cn$lab, 
#                                   levels = sort(unique(go_anno_gene_sel_cn$lab))) # same as defualt

go_anno_gene_sel_cn$group.y = factor(go_anno_gene_sel_cn$group.y, 
                                     levels = c("hummingbirds", "parrots", "honeyeaters", "sunbirds"))

# v2: num cnee per term per clade
go_anno_gene_sel_cn2 = go_anno_gene_sel_cn %>% 
  group_by(go_id, group.y) %>% # per term per clade, how many (unique) cnees, (sum of all genes)
  mutate(n_cnees_perclade_perterm = length(unique(id)))

colp = brewer.pal(4, 'Dark2')
# colp2 = c(rep(colp[1], 2), rep(colp[2], 2), colp[3])

N = 4
size_max = ceiling(max(go_anno_gene_sel_cn2$n_cnees_perclade_perterm)/N)*N
x = go_anno_gene_sel_cn2 %>% filter(go_id =="GO:0051179")

p_term2 = go_anno_gene_sel_cn2 %>% ggplot() +
  geom_point(aes(y = lab, x = group.y, color = group.y, size = n_cnees_perclade_perterm)) +
  facet_grid(ont ~ ., scales = 'free', space = 'free') +
  # scale_x_discrete(name = "# of CNEEs in accelerated state") +
  scale_color_manual(values = colp, guide = "none") +
  scale_size_continuous(name = "Number of accelerated CNEEs", range = c(0.2, 4), limits = c(1, size_max)) +
  # guides() +
  theme_m + 
  theme(#plot.margin = unit(rep(4,4), "pt"), # trbl
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1),
    axis.text.y = element_text(hjust = 1),
    legend.position = "bottom"); p_term2

graph_path = paste0(out_path, "GO_diff_convergent_lv_", tmp_set, "_perterm_ncneeperclade_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*1.3, height = Width_HalfCol*1.5, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p_term2
dev.off()
openfile(graph_path)

# v3: num cnee per term per clade, gene lv
go_anno_gene_sel_cn2 = go_anno_gene_sel_cn %>% 
  group_by(go_id, group.y) %>% # per term per clade, how many (unique) cnees, (sum of all genes)
  mutate(n_cnees_perclade_perterm = length(unique(id))) %>% 
  group_by(go_id, group.y, gene) %>% # per term per clade, how many (unique) cnees, per gene
  mutate(n_cnees_perclade_perterm_pergene = length(unique(id)))

colp = brewer.pal(4, 'Dark2')
# colp2 = c(rep(colp[1], 2), rep(colp[2], 2), colp[3])

N = 4
size_max = ceiling(max(go_anno_gene_sel_cn2$n_cnees_perclade_perterm_pergene)/N)*N
x = go_anno_gene_sel_cn2 %>% filter(go_id =="GO:0051179")

p_term3 = go_anno_gene_sel_cn2 %>% 
  # filter(ont == "BP") %>% 
  # filter(n_clades_per_gene >= 2) %>% 
  ggplot() +
  geom_point(aes(y = gene, x = group.y, color = group.y, size = n_cnees_perclade_perterm_pergene)) +
  facet_grid(ont+lab~.  , scales = 'free', space = 'free') +
  # facet_wrap(lab~.  , scales = 'free', ncol = 3) +
  # scale_x_discrete(name = "# of CNEEs in accelerated state") +
  scale_y_discrete(limits = rev) +
  scale_color_manual(values = colp) + # , guide = "none"
  scale_size_continuous(name = "Number of accelerated CNEEs", range = c(0.2, 4), limits = c(1, size_max)) +
  guides(color = guide_legend(nrow = 4)) +
  theme_m + 
  theme(#plot.margin = unit(rep(4,4), "pt"), # trbl
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1),
    axis.text.y = element_text(hjust = 1),
    strip.text.y = element_text(angle = 0, hjust = 0),
    legend.position = "bottom"); p_term3

graph_path = paste0(out_path, "GO_diff_convergent_lv_", tmp_set, "_perterm_ncneeperclade_genelv_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2.5, height = Width_HalfCol*4, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p_term3
dev.off()
openfile(graph_path)

# E. simplified ####
# by rerun enrichment of the real set
# goals:
# 1. to be able to simplified the GO results and using some clusterprofilers' functions
# 2. get gene ratios

simplifiedGO_func = function(s){
  # s = target_files[1]
  print(paste0("current: ", s))
  tmp_name = basename(s)
  tmp_name = gsub("cnee_goperms_(.*)\\.counts\\.bed", "\\1", tmp_name )
  print(paste0("current: ", tmp_name))
  
  # load in background data and clean up NCBI IDs
  bg <- read_tsv(s, col_names = F, show_col_types = FALSE) %>%
    dplyr::rename(chr = X1, start = X2, end = X3, gene = X4, combo = X5, 
                  biotype = X6, accel = X7, total = X8) %>%
    separate(combo, into = c(NA, "ncbi"), sep = "GeneID:", remove = T) %>%
    dplyr::select(-c(chr, start, end, ncbi, total, biotype))
  
  # set up enrichment calculation
  calc_enrich <- function(targetset, background, tmp_ont) {
    sel_custom_go_anno = go_annoC  %>%
      filter(ont == tmp_ont)
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
  
  # filter out observed target set from background data
  target <- bg %>% filter(accel >= 1)
  
  # calculate enrichment, OR, and enrichment score for real target set for each subontology
  bp.real <- calc_enrich(target, bg, tmp_ont = "biological_process")
  mf.real <- calc_enrich(target, bg, "molecular_function")
  cc.real <- calc_enrich(target, bg, "cellular_component")
  
  return(list('bp' = bp.real, 'mf' = mf.real, 'cc' = cc.real))
}

target_files = list.files(project_path, ".counts.bed$", full.names = T); target_files
target_files = target_files[3]
simplifiedGO = lapply(target_files, simplifiedGO_func)

# e.1 4sets ====
d <- GOSemSim::godata('org.Hs.eg.db', ont="BP")
## y = go_anno_gene %>% filter(group == "4sets")
y = go_anno %>% filter(group == "4sets")
# y = tmp_4sets_bp_top20
# y = tmp_4sets_bp_top50

length(unique(y$go_id)) # 105
y2 = y %>% dplyr::select(go_id, pVal_enrich) %>% distinct()

# enricher - works
# x = simplifiedGO[[3]]
x = simplifiedGO[[1]]
x.bp1 = x[[1]]
x.bp1@organism = "Homo sapiens"
x.bp1@ontology = "BP"
x.bp1@keytype = "SYMBOL"
x.bp1@readable <- TRUE 

# keep only the terms based on permutation and replace padj using the permuted p
x.bp_sel = x.bp1 %>% filter(ID %in% y$go_id)
x.bp_sel2 = x.bp_sel@result %>% left_join(y2, by = c("ID" = "go_id")) %>% 
  dplyr::select(-`p.adjust`) %>% 
  dplyr::rename(`p.adjust` = pVal_enrich)
x.bp_sel@result <- x.bp_sel2 # works
# x.bp_sel_tb = x.bp_sel@result

# x.bp_sel_readable <- setReadable(x.bp_sel, OrgDb = org.Gg.eg.db, keyType="ENTREZID") # works


bp <- enrichplot::pairwise_termsim(x.bp_sel, method = "Wang", semData = d)
nrow(bp@result) 
p1 <- barplot(bp, color = "p.adjust", showCategory=20); p1 # showCategory=20, 
str(p1$data, 1)
top20_termsim = p1$data$ID

p2 <- barplot(bp, color = "p.adjust", showCategory=10); p2 # showCategory=20, 
str(p2$data, 1)
top10_termsim = p2$data$ID

p4 <- barplot(bp, color = "p.adjust", showCategory=30); p2 # showCategory=20, 
top30_termsim = p4$data$ID


bp2 <- enrichplot::pairwise_termsim(x.bp_sel)
nrow(bp2@result) 
p3 <- barplot(bp2, color = "p.adjust", showCategory=20); p3 # showCategory=20, 
str(p3$data, 1)
top20_termsim2 = p3$data$ID

cowplot::plot_grid(p1, p3, ncol=2, labels = LETTERS[1:2]) # doesn't make difference (semData = d)

####

p2 = dotplot(bp, showCategory=15); p2

# p3 <- cnetplot(bp, foldChange=geneList); p3

# p4 <- heatplot(bp, showCategory=5); p4

p5 <- emapplot(bp, showCategory = 25, cex_label_category = 0.6); p5

library(DOSE)
data(geneList)
de <- names(geneList)[abs(geneList) > 2]
edo <- enrichDGN(de)
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')

# hsa
d <- GOSemSim::godata('org.Hs.eg.db', ont="BP")
bp <- enrichplot::pairwise_termsim(x.bp_sel, method = "Wang", semData = d)
nrow(bp@result) # 104 | top20: 20
bp2 <- simplify(bp, cutoff=0.6, by="p.adjust", select_fun=min)
nrow(bp2@result) # 0.6:98
setdiff(bp@result$ID, bp2@result$ID)
bp3 <- simplify(bp, cutoff=0.6, by="p.adjust", select_fun=min)
nrow(bp2@result) # 0.6:98

p1 <- barplot(bp, showCategory=10, color = "p.adjust") 
p2 <- barplot(bp2, showCategory=10, color = "p.adjust") 
cowplot::plot_grid(p1, p2, ncol=2, labels = LETTERS[1:2])


p1 <- treeplot(bp, hilight = F); p1
p2 <- treeplot(bp, hclust_method = "average", 
               hilight = F, offset = rel(50)); p2
aplot::plot_list(p1, p2, tag_levels='A')

emapplot(bp, showCategory = 25, cex_label_category = 0.6)

dotplot(bp)

# ggal
dg <- GOSemSim::godata('org.Gg.eg.db', ont="BP")
bpg <- enrichplot::pairwise_termsim(x.bp_sel, method = "Wang", semData = dg)
nrow(bpg@result) # 104
bpg2 <- simplify(bp, cutoff=0.6, by="p.adjust", select_fun=min)
nrow(bpg2@result) # 0.6:98

p1g <- barplot(bpg, showCategory=10, color = "p.adjust") 
p2g <- barplot(bpg2, showCategory=10, color = "p.adjust") 
cowplot::plot_grid(p1g, p2g, ncol=2, labels = LETTERS[1:2])


p1 <- cnetplot(bp, node_label="all", 
               cex_gene = 0.5,
               cex_label_category = 1); p1 




# "Resnik", "Lin", "Rel", "Jiang" , "Wang" and "JC"
d <- GOSemSim::godata('org.Gg.eg.db', ont="BP")
sapply(c("Resnik", "Lin", "Rel", "Jiang" , "Wang", "JC"), function(m){
  # m = "Lin"
  bp <-  enrichplot::pairwise_termsim(x.bp_sel, method = m, semData = d)
  p2 <- treeplot(bp, hclust_method = "average", 
                 hilight = F, offset = rel(50)); p2
  
  graph_path = paste0(out_path, "treeplot_4sets_", m, "_", Sys.Date(), ".pdf")
  pdf(graph_path, width = Width_HalfCol*3, height = Width_HalfCol*2, pointsize = AxisTxFontSizeSize, onefile = TRUE)
  print(p2)
  dev.off()
  openfile(graph_path)
})


d <- GOSemSim::godata('org.Hs.eg.db', ont="BP")
sapply(c("Resnik", "Lin", "Rel", "Jiang" , "Wang", "JC"), function(m){
  # m = "Lin"
  bp <-  enrichplot::pairwise_termsim(x.bp_sel, method = m, semData = d)
  p2 <- treeplot(bp, hclust_method = "average", 
                 hilight = F, offset = rel(50)); p2
  
  graph_path = paste0(out_path, "treeplot_4sets_", m, "_human_", Sys.Date(), ".pdf")
  pdf(graph_path, width = Width_HalfCol*3, height = Width_HalfCol*2, pointsize = AxisTxFontSizeSize, onefile = TRUE)
  print(p2)
  dev.off()
  openfile(graph_path)
})


### compareCluster [doesn't work] ====
data(gcSample)
xx <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)
dotplot(xx)
library(org.Hs.eg.db)

get_genelist = lapply(target_files[7:10], function(s){
  # s = target_files[7]
  # print(paste0("current: ", s))
  tmp_name = basename(s)
  tmp_name = gsub("cnee_goperms_(.*)\\.counts\\.bed", "\\1", tmp_name )
  print(paste0("current: ", tmp_name))
  
  # load in background data and clean up NCBI IDs
  bg <- read_tsv(s, col_names = F, show_col_types = FALSE) %>%
    dplyr::rename(chr = X1, start = X2, end = X3, gene = X4, combo = X5, 
                  biotype = X6, accel = X7, total = X8) %>%
    separate(combo, into = c(NA, "ncbi"), sep = "GeneID:", remove = T) 
  # target <- bg %>% filter(accel >= 1) %>% pull(ncbi)
  target <- bg %>% filter(accel >= 1) %>% pull(gene)
  x = mapIds(org.Hs.eg.db, target, "ENTREZID", "SYMBOL")
  # x = na.omit(x)
  # return(x)
})

xx <- compareCluster(get_genelist, 
                     # fun="enrichKEGG", organism="gga", # no enrichment
                     # fun = "enrichGO", OrgDb='org.Gg.eg.db' # no enrichment
                     fun = "enrichGO", OrgDb='org.Hs.eg.db' # no enrichment
                     # fun = "enrichPathway",
                     # pvalueCutoff=0.05
                     )
dotplot(xx)

####



# 1: cnee_goperms_2wayconvergent.counts.bed
# y = go_gene_list_reorg %>% filter(group == "2wayconvergent")
# length(unique(y$go_id)) # 22

# 2: cnee_goperms_3wayconvergent.counts.bed
# y = go_gene_list_reorg %>% filter(group == "3wayconvergent")

# 6: genelv_conv_4way
y = go_gene_list_reorg %>% filter(group == "genelv_conv_4way")
length(unique(y$go_id)) # 5

# 5: genelv_conv_4way
# y = go_gene_list_reorg %>% filter(group == "genelv_conv_3way")
# length(unique(y$go_id)) # 24

# enricher - works
x = simplifiedGO[[6]]
x.bp1 = x[[1]]
x.bp1@organism = "Homo sapiens"
x.bp1@ontology = "BP"
x.bp1@keytype = "ENTREZID"

x.bp_sel = x.bp1 %>% filter(ID %in% y$go_id)
x.bp_sel_tb = x.bp_sel@result

x.bp_sel_readable <- setReadable(x.bp_sel, OrgDb = org.Gg.eg.db, keyType="ENTREZID") # works

bp <-  enrichplot::pairwise_termsim(x.bp_sel)
nrow(bp@result) # 24
# bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
nrow(bp2@result) # 22
setdiff(bp@result$ID, bp2@result$ID)

p1 <- barplot(bp, showCategory=25) 
p2 <- barplot(bp2, showCategory=25) 
cowplot::plot_grid(p1, p2, ncol=2, labels = LETTERS[1:2])


# enricherGO
x = simplifiedGO2[[1]]
x.bp2 = x[[1]]

x.bp_sel_readable <- setReadable(x.bp2, OrgDb = org.Gg.eg.db, keyType="ENTREZID") # works

bp <-  enrichplot::pairwise_termsim(x.bp_sel)

bp <-  enrichplot::pairwise_termsim(x.bp)
bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
p1 <- emapplot(bp)
p2 <- emapplot(bp2)
cowplot::plot_grid(p1, p2, ncol=2, labels = LETTERS[1:2])


