# https://github.com/sjswuitchik/duck_comp_gen/blob/master/03_cnee_analyses/05c_go_perms.R

#### run GO analysis ####
library(org.Gg.eg.db) # BiocManager::install("org.Gg.eg.db")
library(tidyverse)
library(clusterProfiler)
library(RColorBrewer)

FontSize = 2.5
AxisTxFontSizeSize = 8
AxisTxFontSizeSize_s = 6
AxisTitleFontSizeSize = 10

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
                axis.line = element_line(color = "black", size = 0.5),
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

project_path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis_2021/05_perm_phyloacc/1_prep_new/postPhyloAcc/go_perms/")
target_files = list.files(project_path, ".counts.bed$", full.names = T); target_files

# do perm ====
do_perm = lapply(target_files, function(s){
  # s = target_files[4]
  print(paste0("current: ", s))
  tmp_name = basename(s)
  tmp_name = gsub("cnee_goperms_(.*)\\.counts\\.bed", "\\1", tmp_name )
  print(paste0("current: ", tmp_name))

  # load in background data and clean up NCBI IDs
  bg <- read_tsv(s, col_names = F) %>%
    rename(chr = X1, start = X2, end = X3, gene = X4, combo = X5, biotype = X6, accel = X7, total = X8) %>%
    separate(combo, into = c(NA, "pass1"), sep = "GeneID:", remove = T) %>%
    separate(pass1, into = c("ncbi", NA), sep = "miRBase:", remove = T) %>%
    select(-c(chr, start, end, gene, total, biotype))
  
  # set up enrichment calculation
  calc_enrich <- function(targetset, background, ont) { 
    enrichGO(targetset$ncbi,'org.Gg.eg.db',
             pvalueCutoff=1.5,
             qvalueCutoff = 1.5, minGSSize=5, maxGSSize=2000, 
             pAdjustMethod="BH",
             universe=background$ncbi,
             keyType="ENTREZID",
             ont=ont) 
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
      select(ncbi, .data[[col_name]])
    bp <- calc_enrich(perm.target, bg, ont = "BP")
    mf <- calc_enrich(perm.target, bg, ont = "MF")
    cc <- calc_enrich(perm.target, bg, ont = "CC")
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
    bp.perms[[i]] <- bp.clean
    mf.perms[[i]] <- mf.clean
    cc.perms[[i]] <- cc.clean
    
  }
  
  # bind permutation results
  bp.perms.clean <- bind_rows(list(bp.perms), .id = "perm")
  mf.perms.clean <- bind_rows(list(mf.perms), .id = "perm")
  cc.perms.clean <- bind_rows(list(cc.perms), .id = "perm")
  
  # write out perms so you don't have to re-run when your session inevitably crashes
  path = paste0(project_path, tmp_name, "_bp.perms.clean.tsv")
  write_tsv(bp.perms.clean, path, col_names = T, na = "")
  path = paste0(project_path, tmp_name, "_mf.perms.clean.tsv")
  write_tsv(mf.perms.clean, path, col_names = T, na = "")
  path = paste0(project_path, tmp_name, "_cc.perms.clean.tsv")
  write_tsv(cc.perms.clean, path, col_names = T, na = "")
  
  # filter out observed target set from background data
  target <- bg %>% filter(accel >= 1)
  
  # calculate enrichment, OR, and enrichment score for real target set for each subontology
  bp.real <- calc_enrich(target, bg, "BP")
  mf.real <- calc_enrich(target, bg, "MF")
  cc.real <- calc_enrich(target, bg, "CC")
  
  bp.real.clean <- bp.real@result %>% 
    separate(GeneRatio, into = c("target_in", "target_total"), remove = F) %>%
    separate(BgRatio, into = c("bg_in", "bg_total")) %>%
    mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
           bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
           OR = target_frac/bg_frac,
           enrich = log2(target_frac/bg_frac),
           newpval = ifelse(is.na(pvalue), 1, pvalue),
           logp = -log10(newpval)) %>%
    dplyr::select(ID, logp, target_frac, target_total, bg_frac, enrich) %>%
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
    dplyr::select(ID, logp, target_frac, target_total, bg_frac, enrich) %>%
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
    dplyr::select(ID, logp, target_frac, bg_frac, enrich) %>%
    arrange(ID)
  
  # merge observed results with permutation results
  bp.merge <- bp.real.clean %>% 
    left_join(., bp.perms.clean, by = c("ID" = "ID"), suffix = c(".real", ".perm")) %>% 
    distinct() %>%
    arrange(ID) %>%
    select(-perm)
  mf.merge <- mf.real.clean %>% 
    left_join(., mf.perms.clean, by = c("ID" = "ID"), suffix = c(".real", ".perm")) %>% 
    distinct() %>%
    arrange(ID) %>%
    select(-perm)
  cc.merge <- cc.real.clean %>% 
    left_join(., cc.perms.clean, by = c("ID" = "ID"), suffix = c(".real", ".perm")) %>% 
    distinct() %>%
    arrange(ID) %>%
    select(-perm)
  
  # calculate empirical pvalues (number of instances where real data is greater than or equal to permutation)
  bp.p <- bp.merge %>%
    group_by(ID) %>%
    mutate(gt_target = target_frac.perm <= target_frac.real,
           gt_enrich = enrich.perm <= enrich.real)
  bp.cols <- sapply(bp.p[,10:11], as.numeric)
  bp.pval <- bp.p %>%
    group_by(ID) %>%
    mutate(sum_target = sum(gt_target) + 1,
           pVal_target = sum_target/1001,
           sum_enrich = sum(gt_enrich) + 1,
           pVal_enrich = sum_enrich/1001) %>%
    select(ID, pVal_target, pVal_enrich) %>%
    distinct()
  
  mf.p <- mf.merge %>%
    group_by(ID) %>%
    mutate(gt_target = target_frac.perm <= target_frac.real,
           gt_enrich = enrich.perm <= enrich.real)
  mf.cols <- sapply(mf.p[,10:11], as.numeric)
  mf.pval <- mf.p %>%
    group_by(ID) %>%
    mutate(sum_target = sum(gt_target) + 1,
           pVal_target = sum_target/1001,
           sum_enrich = sum(gt_enrich) + 1,
           pVal_enrich = sum_enrich/1001) %>%
    select(ID, pVal_target, pVal_enrich) %>%
    distinct()
  
  cc.p <- cc.merge %>%
    group_by(ID) %>%
    mutate(gt_target = target_frac.perm <= target_frac.real,
           gt_enrich = enrich.perm <= enrich.real)
  cc.cols <- sapply(cc.p[,10:11], as.numeric)
  cc.pval <- cc.p %>%
    group_by(ID) %>%
    mutate(sum_target = sum(gt_target) + 1,
           pVal_target = sum_target/1001,
           sum_enrich = sum(gt_enrich) + 1,
           pVal_enrich = sum_enrich/1001) %>%
    select(ID, pVal_target, pVal_enrich) %>%
    distinct()
  
  # write out GO terms
  path = paste0(project_path, tmp_name, "_bp_GOterms.tsv")
  write_tsv(bp.pval, path, col_names = T, na = "")
  path = paste0(project_path, tmp_name, "_mf_GOterms.tsv")
  write_tsv(mf.pval, path, col_names = T, na = "")
  path = paste0(project_path, tmp_name, "_cc_GOterms.tsv")
  write_tsv(cc.pval, path, col_names = T, na = "")
  
  
  # write out sig GO terms
  bp.pval_sig = filter(bp.pval, pVal_enrich <= 0.05) %>% 
    dplyr::select(-c(pVal_target)) %>% mutate(subontology = "BP") 
  path = paste0(project_path, tmp_name, "_bp_GOterms_sig.tsv")
  write_tsv(bp.pval_sig, path, col_names = T, na = "")
  
  mf.pval_sig = filter(mf.pval, pVal_enrich <= 0.05) %>% 
    dplyr::select(-c(pVal_target)) %>% 
    mutate(subontology = "MF") 
  path = paste0(project_path, tmp_name, "_mf_GOterms_sig.tsv")
  write_tsv(mf.pval_sig, path, col_names = T, na = "")
  
  cc.pval_sig = filter(cc.pval, pVal_enrich <= 0.05) %>% 
    dplyr::select(-c(pVal_target)) %>% 
    mutate(subontology = "CC")
  path = paste0(project_path, tmp_name, "_cc_GOterms_sig.tsv")
  write_tsv(cc.pval_sig, path, col_names = T, na = "")
  
  sig_comb = bp.pval_sig %>% bind_rows(mf.pval_sig) %>% bind_rows(cc.pval_sig) 
  path = paste0(project_path, tmp_name, "_combined_GOterms_sig.tsv")
  write_tsv(sig_comb, path, col_names = T, na = "")
  
  return(sig_comb)
  # bp.perms.clean, mf.perms.clean, cc.perms.clean
  # bp.pval, mf.pval, cc.pval,
  # bp.pval_sig, mf.pval_sig, cc.pval_sig,
  # sig_comb
})
  

# analysis ####
library(biomaRt) # installed dev version with BiocManager::install('grimbough/biomaRt')

sig_files = list.files(project_path, "_combined_GOterms_sig.tsv", full.names = T); sig_files

out2 = lapply(sig_files, function(s){
  # s = sig_files[2]
  tmp_name = basename(s)
  sp = gsub("(\\w+)_combined_GOterms_sig.tsv", "\\1", tmp_name)
  print(sp)
  goList <- read_tsv(s, col_names = T, col_types = c('c', 'n', 'c')) %>% 
    # dplyr::select(ID) %>% 
    bind_cols('group' = sp)
}) 
out3 = out2 %>% bind_rows() %>%
  rename(goID = ID)

#### annotate GO IDs with functional annotation
goList = out3$goID
mart <- useMart(biomart = 'ensembl', dataset = 'ggallus_gene_ensembl')
martList <- getBM(attributes = c("external_gene_name", "go_id", "name_1006"), values = goList, bmHeader = T, mart = mart)

#### pull name, not gene 
collapse <- martList %>% 
  dplyr::rename(goID= `GO term accession`, goTerm = `GO term name`) %>%
  group_by(goID) %>%
  summarise(goTerm = paste(sort(unique(goTerm)), collapse = ", ")) 

clean_mart = collapse %>%
  # remove rows without goID
  filter(! (is.na(goID) | goID == "" )) %>%
  # remove rows without goTerm
  filter(! (is.na(goTerm) | goTerm == "" ))

go_anno = out3 %>% left_join(clean_mart, by = "goID") %>% 
  dplyr::select(goID, group, subontology, pVal_enrich, goTerm) %>%
  arrange(group, subontology) %>%
  mutate(goTerm = ifelse(grepl("GO:0048518", goID), "positive regulation of biological process", goTerm),
         goTerm = ifelse(grepl("GO:0006928", goID), "movement of cell or subcellular component", goTerm),
         goTerm = ifelse(grepl("GO:0022610", goID), "obsolete biological adhesion", goTerm),
         goTerm = ifelse(grepl("GO:0032501", goID), "multicellular organismal process", goTerm),
         goTerm = ifelse(grepl("GO:0051641", goID), "cellular localization", goTerm),
         goTerm = ifelse(grepl("GO:0051674", goID), "localization of cell", goTerm),
         goTerm = ifelse(grepl("GO:0120035", goID), "regulation of plasma membrane bounded cell projection organization", goTerm),
         goTerm = ifelse(grepl("GO:0120036", goID), "plasma membrane bounded cell projection organization", goTerm),
         goTerm = ifelse(grepl("GO:0065007", goID), "biological regulation", goTerm),
         goTerm = ifelse(grepl("GO:0005622", goID), "intracellular anatomical structure", goTerm),
         goTerm = ifelse(grepl("GO:0043228", goID), "non-membrane-bounded organelle", goTerm)
         )

go_anno%>% filter(is.na(goTerm)) %>% select(goID) %>% distinct()

path = paste0(project_path, "Final_GO_combined_sig_mart_anno.tsv")
# write_tsv(go_anno, path, col_names = T, na = "")

go_anno_tmp = go_anno %>%
  mutate(lab = paste0(goID, ":", goTerm), 
         `-log10P` = -log10(pVal_enrich))
go_anno_tmp$group = factor(go_anno_tmp$group, levels = c("hummingbirds", "parrots", 
                                                         "honeyeaters_pardalote", "sunbirds_flowerpecker"))
go_anno_tmp$subontology = factor(go_anno_tmp$subontology, levels = c("BP", "MF", "CC"))

pal = brewer.pal(n = 4, name = "Dark2")

pal1 = pal[1]
p1 = ggplot(go_anno_tmp %>% filter(group == "hummingbirds") ) +
  geom_col(aes(x = `-log10P`, y = lab, fill = group), alpha = 0.7) +
  scale_x_continuous(expand = c(0,0,0.05,0)) +
  facet_grid(facets = group+subontology~., space = 'free', scale = 'free') + 
  scale_fill_manual(values = pal1) +
  theme_m + 
  theme(axis.title.y = element_blank(), 
        legend.position = "bottom"); p1
graph_path = paste0(project_path, "hummingbirds_GO_", Sys.Date(), ".pdf")
pdf(graph_path, width = 7.7, height = 5.85, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p1
dev.off()


pal2 = pal[-1]
p2 = ggplot(go_anno_tmp %>% filter(group != "hummingbirds") ) +
  geom_col(aes(x = `-log10P`, y = lab, fill = group), alpha = 0.7) +
  scale_x_continuous(expand = c(0,0,0.05,0)) +
  facet_grid(facets = group+subontology~., space = 'free', scale = 'free') + 
  scale_fill_manual(values = pal2) +
  theme_m + 
  theme(axis.title.y = element_blank(), 
        legend.position = "bottom"); p2
graph_path = paste0(project_path, "except_GO_", Sys.Date(), ".pdf")
pdf(graph_path, width = 7.7, height = 5.85, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p2
dev.off()


