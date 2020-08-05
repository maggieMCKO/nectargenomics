#revised Nov 2018 for paper revisions
library(tidyverse)
library(cluster)
library(clusterProfiler)
library(org.Gg.eg.db)

path_to_data <- paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/2_GO/1_enrichment/")

# 1. ori script ====
transform_real <- function(DF) {
  DF %>% 
    separate(GeneRatio, into=c("target_in", "target_total")) %>% 
    separate(BgRatio, into=c("bg_in", "bg_total")) %>%
    mutate(newpval = ifelse(is.na(pvalue), 1, pvalue), 
           logp = -log10(newpval),
           target_frac = as.numeric(target_in)/as.numeric(target_total), 
           bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
           enrich = log2(target_frac/bg_frac)) %>%
    dplyr::select(version, set, ID, Description, logp, target_frac, bg_frac, enrich, geneID, Count) %>% 
    filter(target_frac < 1, bg_frac < 1) %>%
    arrange(ID)
}

# the real outputs of different runs are all the same, just use the first run
orig_bp <- read_tsv(paste0(path_to_data, "/goperms/original_GO_galgal6_run1_BP_real.tsv")) %>% transform_real()
orig_mf <- read_tsv(paste0(path_to_data, "/goperms/original_GO_galgal6_run1_MF_real.tsv")) %>% transform_real()

#read perms
orig_bp_perm<-readRDS(paste0(path_to_data, "/goperms/original_galgal6_BP.robj"))
orig_mf_perm<-readRDS(paste0(path_to_data, "/goperms/original_galgal6_MF.robj"))

orig_bp_merge <- full_join(orig_bp, orig_bp_perm, by=c("version" = "version", "set" = "set", "ID" = "ID")) 
orig_mf_merge <- full_join(orig_mf, orig_mf_perm, by=c("version" = "version", "set" = "set", "ID" = "ID"))

orig_bp_perm_raw = read_csv(paste0(path_to_data, "/goperms/original_galgal6_BP_raw_2020-07-28.csv")) 
orig_mf_perm_raw = read_csv(paste0(path_to_data, "/goperms/original_galgal6_MF_raw_2020-07-28.csv")) 


#compute P-values

orig_bp_merge <-
  orig_bp_merge %>%
  rowwise %>%
  mutate(pval_frac = max(1-ecdf_frac(target_frac), 0.0002), # if smaller than 0.0002, take 0.0002, don't want 0
         pval_logp = max(1-ecdf_logp(logp), 0.0002), # 1-orig_bp_merge$ecdf_frac[[1]](orig_bp_merge$target_frac[1])
         pval_enrich = max(1-ecdf_enrich(enrich), 0.0002),
         enrichX = ecdf_enrich(enrich)) %>%
  ungroup %>% group_by(version, set) %>%
  mutate(qval_logp = p.adjust(pval_logp, "BH"),
         qval_frac = p.adjust(pval_frac, "BH"),
         qval_enrich = p.adjust(pval_enrich, "BH"))


orig_mf_merge <-
  orig_mf_merge %>% 
  rowwise %>% 
  mutate(pval_frac = max(1-ecdf_frac(target_frac), 0.0002), 
         pval_logp = max(1-ecdf_logp(logp), 0.0002), 
         pval_enrich = max(1-ecdf_enrich(enrich), 0.0002)) %>% 
  ungroup %>% group_by(version, set) %>% 
  mutate(qval_logp = p.adjust(pval_logp, "BH"),
         qval_frac = p.adjust(pval_frac, "BH"),
         qval_enrich = p.adjust(pval_enrich, "BH"))


#analysis

mf = orig_mf_merge %>% filter(version == "basic") %>% filter(qval_frac < 0.25) %>% 
  mutate(exp_frac = map(ecdf_frac, summary) %>% map_dbl(., 3)) %>%  # take median (the 3rd element)
  # select(set, ID, target_frac, exp_frac, bg_frac, enrich, pval_frac, qval_frac) %>% 
  dplyr::select(set, ID, Description, target_frac, exp_frac, bg_frac, enrich, pval_frac, qval_frac, geneID, Count) %>%
  arrange(qval_frac) #%>% mutate(term = go2term(ID)) 
# mf_term =  go2term(mf$ID)
# mf = mf %>% full_join(mf_term, by = c("ID" = "go_id"))  

path = paste0(path_to_data, "/GOPERM_orig_mf_results_", Sys.Date(), ".tsv")  
# write_tsv(mf, path = path)

bp = orig_bp_merge %>%filter(version == "basic") %>% filter(qval_frac < 0.25) %>% 
  mutate(exp_frac = map(ecdf_frac, summary) %>% map_dbl(., 3)) %>%
  # select(set, ID, target_frac, exp_frac, bg_frac, enrich, pval_frac, qval_frac) %>% 
  dplyr::select(set, ID, Description, target_frac, exp_frac, bg_frac, enrich, pval_frac, qval_frac, geneID, Count) %>%
  arrange(qval_frac) 
# bp_term =  go2term(bp$ID)
# bp = bp %>% full_join(bp_term, by = c("ID" = "go_id"))  

path = paste0(path_to_data, "/GOPERM_orig_bp_results_", Sys.Date(), ".tsv")
# write_tsv(bp, path)

## 2. annotate back to go output ====
library(biomaRt) # bioconductor

path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/2_GO/0_prepare/galGal6_final_merged_CNEEs_cloest_genes.bed")
# gene_gg = read_tsv(paste0("../04_wga/03_ce_annotation/cnees.", args[1], ".annotation"), col_names = c("cnee", "gene"))
gene_gg = read_tsv(path, col_names = F)
gene_gg = gene_gg[, c(4, 8)]
names(gene_gg) = c("cnee", "gene")
gene_gg$gene = gsub("gene-", "", gene_gg$gene)

## get gene symbol ====

split_gene = function(input_df){
  # input_df = bp
  df = lapply(1:nrow(input_df), function(n_row){
    tmp_row = input_df[n_row, ];
    out = unlist(strsplit(tmp_row$geneID, "/")); 
    genesymb = bitr(out, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Gg.eg.db);
    term =  tmp_row$Description ;
    goid = tmp_row$ID;
    genesymb = genesymb %>% bind_cols(Description = term, ID = goid)}) 
  }

go_gene_list = lapply(list(bp, mf), split_gene) %>% bind_rows(.id = "gotype")

# use biomaRt to get gene description
# biomaRt https://www.bioconductor.org/packages/devel/bioc/vignettes/biomaRt/inst/doc/biomaRt.html
listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
head(datasets) 
dim(datasets) # 203 3
datasets$dataset[grep("gal", datasets$dataset, ignore.case = T)]

ensembl_gg = useDataset("ggallus_gene_ensembl", mart=ensembl)

filters = listFilters(ensembl_gg) # query # from
filters$name[grep('entr', filters$name)] # entrezgene_id

attributes = listAttributes(ensembl_gg) # to
attributes$name[grep('name', attributes$name)] # wikigene_name external_gene_name
attributes$name[grep('ggallus', attributes$name)] # 

x =  getBM(attributes=c('entrezgene_id', 'description'), 
      filters = 'entrezgene_id', 
      values = go_gene_list$ENTREZID, 
      mart = ensembl_gg)
x$entrezgene_id = as.character(x$entrezgene_id)

go_gene_list2 = go_gene_list %>% dplyr::full_join(x, by = c("ENTREZID" = "entrezgene_id")) %>%  mutate(gotype = ifelse(grepl(1, gotype), "bp", "mf")) %>%
  mutate(geneName = gsub(" \\[Source:NCBI.*", "", description)) %>%
  dplyr::select(gotype, ENTREZID, SYMBOL, geneName, ID, Description) 

go_gene_list_uniq = go_gene_list2 %>% group_by(SYMBOL) %>%
  mutate(Terms = paste(Description, collapse = "; ")) %>% distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::select(-gotype, -Description) # 58

path = paste0(path_to_data, "/GOPERM_orig_combined_results_", Sys.Date(), ".tsv")
# write_tsv(go_gene_list2, path)

path = paste0(path_to_data, "/GOPERM_orig_combined_results_uniq_", Sys.Date(), ".tsv")
# write_tsv(go_gene_list_uniq, path)

# which cnees were associated 
go_gene_list_uniq_cnee = go_gene_list_uniq %>% 
  left_join(gene_gg, by = c("SYMBOL" = "gene")) %>% group_by(SYMBOL) %>% 
  mutate(count_cnee = n()) 

path = paste0(path_to_data, "/GOPERM_orig_combined_results_uniq_cnee_long_", Sys.Date(), ".tsv")
# write_tsv(go_gene_list_uniq_cnee, path)

# collapse cnees
go_gene_list_uniq_cnee2 = go_gene_list_uniq_cnee %>% mutate(cnees = paste(cnee, collapse = "; ")) %>% 
  distinct(SYMBOL, .keep_all = TRUE) %>% dplyr::select(-cnee)

path = paste0(path_to_data, "/GOPERM_orig_combined_results_uniq_cnee_", Sys.Date(), ".tsv")
# write_tsv(go_gene_list_uniq_cnee2, path)

# top10 GO ====

go = bp %>% bind_rows(mf, .id = 'gotype') %>% mutate(gotype = ifelse(grepl(1, gotype), "bp", "mf"))


top10 = go %>% group_by(gotype) %>% 
  mutate(log10_p = -log10(pval_frac)) %>%
  slice_max(order_by = log10_p, n = 10)

go_gene_list_sub = go_gene_list2 %>% filter(ID %in% top10$ID)

go_gene_list_sub_uniq = go_gene_list_sub %>% group_by(SYMBOL) %>%
  mutate(Terms = paste(Description, collapse = "; ")) %>% 
  distinct(SYMBOL, .keep_all = TRUE) %>% dplyr::select(-gotype, -Description) # 31

path = paste0(path_to_data, "/GOPERM_orig_combined_results_top10_", Sys.Date(), ".tsv")
# write_tsv(go_gene_list_sub, path)

path = paste0(path_to_data, "/GOPERM_orig_combined_results_top10_uniq_", Sys.Date(), ".tsv")
# write_tsv(go_gene_list_sub_uniq, path)

# which cnees were associated 
go_gene_list_sub_uniq_cnee = go_gene_list_sub_uniq %>% 
  left_join(gene_gg, by = c("SYMBOL" = "gene")) %>% group_by(SYMBOL) %>% 
  mutate(count_cnee = n()) 

path = paste0(path_to_data, "/GOPERM_orig_combined_results_top10_uniq_cnee_long_", Sys.Date(), ".tsv")
# write_tsv(go_gene_list_sub_uniq_cnee, path)

go_gene_list_sub_uniq_cnee2 = go_gene_list_sub_uniq_cnee %>% 
  mutate(cnees = paste(cnee, collapse = "; ")) %>% 
  distinct(SYMBOL, .keep_all = TRUE) %>% dplyr::select(-cnee)

path = paste0(path_to_data, "/GOPERM_orig_combined_results_top10_uniq_cnee_", Sys.Date(), ".tsv")
# write_tsv(go_gene_list_sub_uniq_cnee2, path)

## compare with eye-balling cnees ====
library(readxl)
path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/1_runphyloacc/1_run_phyloacc/plots2/Acc_cnees.xlsx")
phyloacc_cnee = read_excel(path)

# go_gene_list_uniq = 

x = sapply(phyloacc_cnee$gene, function(s){
  grep(s, go_gene_list_uniq$SYMBOL)
})
x = unlist(x) # CENPW

# 3. plot ====
library(grid)
library(gridExtra)
FontSize = 3
AxisTxFontSizeSize = 8
AxisTxFontSizeSize_s = 5
AxisTitleFontSizeSize = 10
Factor_mmtoin = 0.0393701
Width_HalfCol = 85*Factor_mmtoin 

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

theme_m = 
  theme(
    plot.background = element_blank(),  #element_rect(fill = "transparent"), # defualt is white
    plot.title = element_text(hjust = 0.5, size = AxisTitleFontSizeSize),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(), #element_rect(color = "black", fill= NA, size = 0.1), # 
    panel.spacing.y = unit(0.2, "lines"),
    axis.line = element_line(color = "black", size = 0.1),
    axis.title.x = element_text(colour = "black", size = AxisTitleFontSizeSize),
    axis.title.y = element_text(colour = "black", size = AxisTitleFontSizeSize), # , hjust = 1, vjust = 0.5
    axis.text.x = element_text(colour = "black", size = AxisTxFontSizeSize),
    axis.text.y = element_text(colour = "black", size = AxisTxFontSizeSize), # , hjust = 1
    axis.ticks = element_line(colour = "black", size = 0.05),
    strip.text = element_text(colour = "black", size = AxisTxFontSizeSize, 
                              margin = margin(2,0,2,0, "pt")),
    legend.key = element_blank(),
    legend.text = element_text(colour = "black", size = AxisTxFontSizeSize),
    legend.title = element_text(colour = "black", size = AxisTitleFontSizeSize),
    legend.position = "none")

project_path = path_to_data

# 3.1 plot topGO ====
plot_topGO = function(input_df, fill_col, title){
  # input_df = bp
  # fill_col = 'red'
  input_df = input_df %>% mutate(log10_p = -log10(pval_frac)) %>%
    slice_max(order_by = log10_p, n = 10)
  # arrange(log10_p)
  input_df$Description <- factor(input_df$Description) %>%
    fct_reorder(input_df$log10_p)
  p = ggplot(input_df, aes(x = log10_p, y = Description)) +
    geom_col(fill = alpha(fill_col, 0.35)) +
    geom_text(aes(label = Count, x = log10_p*1.1), size = FontSize) +
    scale_x_continuous(name = "-log10(pval)") +
    ggtitle(title) +
    theme_m + theme(axis.title.y = element_blank()); p
}

plot_topGO(bp, 'blue')

p_t = mapply(plot_topGO, input_df = list(bp, mf), fill_col= c('red', 'blue'), title = c('Biological process', 'Molecular function'), SIMPLIFY = FALSE)
p = arrangeGrob(grobs = p_t, widths = c(4, 5))

graph_path = paste0(project_path, "Top10GOterms_", Sys.Date(), ".pdf")       # all data
pdf(graph_path, width = Width_HalfCol*3, height = Width_HalfCol*0.8, pointsize = AxisTxFontSizeSize, onefile = TRUE)
# png(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*0.625, unit = "in", res = 600)
grid.draw(p)
dev.off()
openfile(graph_path)

# 3.2 plot topGENE by associated CNEE ====
plot_topGene = function(input_df){
  # input_df = go_gene_list_sub_uniq_cnee2
  fill_col = 'red'
  input_df = input_df %>% slice_max(order_by = count_cnee, n = 10)
  # arrange(log10_p)
  input_df$SYMBOL <- factor(input_df$SYMBOL) %>%
    fct_reorder(input_df$count_cnee)
  
  p = ggplot(input_df, aes(x = count_cnee, y = SYMBOL)) +
    geom_col(fill = alpha(fill_col, 0.35)) +
    # geom_text(aes(label = count_cnee, x = count_cnee*1.1), size = FontSize) +
    scale_x_continuous(name = "# associated CNEEs") +
    # ggtitle(title) +
    theme_m + theme(axis.title.y = element_blank()); p
}

p_t = lapply(list(go_gene_list_uniq_cnee2, go_gene_list_sub_uniq_cnee2), plot_topGene)
p = arrangeGrob(grobs = p_t, widths = c(4, 5))

graph_path = paste0(project_path, "NumClosestCnee_byGene_", Sys.Date(), ".pdf")       # all data
pdf(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*2, pointsize = AxisTxFontSizeSize, onefile = TRUE)
# png(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*0.625, unit = "in", res = 600)
grid.draw(p)
dev.off()
openfile(graph_path)

## 3.3 perm vs real [top10] ====

plot_enrich = function(input_perm_raw, tmp_gotype){
  
  top10_gptype = top10 %>% filter(gotype == tmp_gotype)  %>% 
    mutate(des = paste0(ID, "\n", Description))
  top10_gptype_sel = top10_gptype %>% dplyr::select(ID, Description, qval_frac, des)
  
  top10_raw_gotype = input_perm_raw  %>% mutate(p_ori = 1/10^logp.perm, qval = p.adjust(p_ori, "BH")) %>%
    # filter(ID %in% top10$ID) %>% 
    right_join(top10_gptype_sel, by = 'ID')
  # dim(top10_raw_gotype) # bp: 34301    10
  # range(top10_raw_gotype$enrich) # bp: -2.086565  3.640936
  # range(top10_raw_gotype$qval)   # bp: 0.0004305326 0.9952342510
  
  top10_raw_gotype$ID = factor( top10_raw_gotype$ID) %>% fct_reorder(top10_raw_gotype$qval_frac)
  
  p = ggplot() +
    geom_density(data = top10_raw_gotype, aes(enrich)) +
    geom_vline(data = top10_gptype, aes(xintercept = enrich), color = 'red') +
    facet_grid(.~des) + theme_m + theme(strip.text = element_text(size = AxisTxFontSizeSize_s)); p
}

p_out = mapply(plot_enrich, input_perm_raw = list(orig_bp_perm_raw, orig_mf_perm_raw), tmp_gotype = c('bp', 'mf'), SIMPLIFY = F)
p = arrangeGrob(grobs = p_out, nrow = 2)

dev.off()
grid.draw(p)

graph_path = paste0(project_path, "PermEnrich_", Sys.Date(), ".pdf")       # all data
pdf(graph_path, width = Width_HalfCol*5, height = Width_HalfCol*1, pointsize = AxisTxFontSizeSize, onefile = TRUE)
# png(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*0.625, unit = "in", res = 600)
grid.draw(p)
dev.off()
openfile(graph_path)

ggplot() +
  geom_density(data =  x, aes(x)) 
  # geom_vline(data = top10_bp, aes(xintercept = enrich), color = 'red') +
  # facet_grid(.~ID) + theme_m
