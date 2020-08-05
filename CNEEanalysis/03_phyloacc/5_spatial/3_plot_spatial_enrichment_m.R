#plot spatial enrichment - clean up to do right

library(tidyverse)
library(qqman)


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

# setwd("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/")
project_path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/5_spatial/")

path = paste0(project_path, "original_spatial_results.tsv")
orig = read_tsv(path, guess_max = 10000) %>% filter(version=="basic")

# orig_gm <- read_tsv("original_spatial_results.tsv.gz", guess_max = 10000) %>% filter(version=="gain")
# orig_gl <- read_tsv("original_spatial_results.tsv.gz", guess_max = 10000) %>% filter(version=="gain_gap")
# ext_gm <- read_tsv("extended_spatial_results.tsv.gz", guess_max = 10000) %>% filter(version=="gain")
# ext_gl <- read_tsv("extended_spatial_results.tsv.gz", guess_max = 10000) %>% filter(version=="gain_gap")
# red_gm <- read_tsv("reduced_spatial_results.tsv.gz", guess_max = 100000) %>% filter(version=="gain")
# red_gl <- read_tsv("reduced_spatial_results.tsv.gz", guess_max = 100000) %>% filter(version=="gain_gap")

#PLOTTING ORIGINAL##
# gwide <- ext_gm %>% filter(rar_qval < 0.01, window_size == "1000kb_100kb_slide") %>% summarize(maxpval=max(rar_pval)) %>% pull(maxpval)
# ss <- ext_gm %>% filter(rar_qval < 0.05, window_size == "1000kb_100kb_slide") %>% summarize(maxpval = max(rar_pval)) %>% pull(maxpval)
gwide <- orig %>% filter(rar_qval < 0.01, window_size == "1000kb_100kb_slide") %>% summarize(maxpval=max(rar_pval)) %>% pull(maxpval)
ss <- orig %>% filter(rar_qval < 0.05, window_size == "1000kb_100kb_slide") %>% summarize(maxpval = max(rar_pval)) %>% pull(maxpval)

plot_func = function(window_size_tmp){
  
  gwide <- orig %>% filter(rar_qval < 0.01, window_size == window_size_tmp) %>% summarize(maxpval=max(rar_pval)) %>% pull(maxpval)
  ss <- orig %>% filter(rar_qval < 0.05, window_size == window_size_tmp) %>% summarize(maxpval = max(rar_pval)) %>% pull(maxpval)
  
  #convert to format for qqman
  unique(orig$chr)
  se_manh <- orig %>% filter(window_size == window_size_tmp) %>% select(chr, start, end, rar_pval, window) %>%
    mutate(chr = sub("chr", "", chr)) %>% 
    # mutate(CHR = as.numeric(sub("Z", 29, chr))) %>%
    mutate(CHR = gsub("Z", 29, chr)) %>%
    mutate(CHR = gsub("W", 30, CHR)) %>%
    mutate(BP = start+((end+1-start)/2)) %>% 
    dplyr::rename(SNP = window) %>%
    dplyr::rename(P = rar_pval) %>% select(SNP, CHR, BP, P)
  se_manh$CHR = as.numeric(se_manh$CHR)
  unique(se_manh$CHR)
  
  #write plot to pdf
  
  # pdf("fig3c.pdf")
  graph_path = paste0(project_path, "manhattan_", window_size_tmp, ".pdf")
  pdf(graph_path, width = Width_HalfCol*3, height = Width_HalfCol*1.5, pointsize = AxisTxFontSizeSize, onefile = TRUE)
  manhattan(se_manh, chrlabs=c(1:15, 17:28, "Z", "W", 33), suggestiveline=-log10(ss), genomewideline = -log10(gwide), col=c("blue4", "orange3"), type="l") 
  dev.off()
  openfile(graph_path)
  
}

tmp = unique(orig$window_size)
lapply(tmp, plot_func)

#windows with an excess of convergent rars - manually annotate in figure
x = orig %>% filter(window_size == "1000kb_100kb_slide") %>% 
  mutate(RAR = rar_pval < gwide, 
         filter_non_target_5sp = filter_non_target_5sp_pval < gwide, 
         filter_non_target_humm2anc = filter_non_target_humm2anc_pval < gwide, 
         filter_non_target_humm2 = filter_non_target_humm2_pval < gwide, 
         filter_non_target_humm2anc_triMol = filter_non_target_humm2anc_triMol_pval < gwide, 
         filter_non_target_humm2anc_honeyanc = filter_non_target_humm2anc_honeyanc_pval < gwide, 
         filter_non_target_honey2anc_triMol = filter_non_target_honey2anc_triMol_pval < gwide, 
         logp = -log10(rar_pval)) %>% 
  # filter(RAR) %>% 
  filter(RAR | filter_non_target_5sp | filter_non_target_humm2anc | filter_non_target_humm2 |
         filter_non_target_humm2anc_triMol | filter_non_target_humm2anc_honeyanc | filter_non_target_honey2anc_triMol) %>%
  # select(chr, start, end, RAR, CRAR, logp) %>%
  select(chr, start, end, logp, RAR,
         filter_non_target_5sp, filter_non_target_humm2anc, filter_non_target_humm2,
         filter_non_target_humm2anc_triMol, filter_non_target_humm2anc_honeyanc, filter_non_target_honey2anc_triMol) 

sum(x$RAR) # 184
sum(x$filter_non_target_5sp) # 0
sum(x$filter_non_target_humm2anc) # 89
sum(x$filter_non_target_humm2) # 124
sum(x$filter_non_target_humm2anc_triMol) # 1
sum(x$filter_non_target_humm2anc_honeyanc) # 0
sum(x$filter_non_target_honey2anc_triMol) # 0

path = paste0(project_path, "gwide_sig_windows.tsv")
# write_tsv(x, path)

# then use bedtools
path = paste0(project_path, "gwide_sig_windows_cloest_genes.bed")
xg = read_tsv(path, col_names = F) # nrow 4473
unique(xg$X1)
xg_m = xg %>% select(X1,X2,X3,X4,X6,X7,X8) %>%
  separate(X8, into = c("useless", 'sym'), sep = "-") %>% 
  select(-useless) %>% 
  filter(X4 > -log10(gwide)) %>%
  arrange(desc(X4), sym) %>%
  group_by(X1, sym) %>% distinct(sym, .keep_all = T)  # nrow 4449 -> 979 (rm group_by X4) -> 916 (X4 > -log10(gwide))
unique(xg_m$X1)

xg_sub = xg_m %>% filter(X1 =='chr3', X4 > 14.76) %>% distinct(sym, .keep_all = T) 
xg_sub = xg_m %>% filter(X1 =='chr5', X4>35.36) %>% distinct(sym, .keep_all = T) 
xg_sub = xg_m %>% filter(X1 =='chr7', X4>19) %>% distinct(sym, .keep_all = T) 
xg_sub = xg_m %>% filter(X1 =='chr14', X4>13) %>% distinct(sym, .keep_all = T) 
xg_sub = xg_m %>% filter(X1 =='chrZ', X4>6) %>% distinct(sym, .keep_all = T) 
cat(xg_sub$sym, sep = "\n")


xg_sub = xg_m %>% filter(X4 > 6) %>% group_by(X1) %>% slice_max(order_by = X4) %>% 
  arrange(desc(X4), X1, sym) 
