#revised Oct 2018 for paper revisions
library(tidyverse)
library(RColorBrewer)

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


# path_to_data <- paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/4_gene/geneperms_first_try/")
# path_to_data <- paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/4_gene/geneperms_2/")
path_to_data <- paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/4_gene/geneperms/")
project_path = path_to_data

#read real data
path = paste0(path_to_data, "original_gene_galgal6_run1_real.tsv")
orig <-read_tsv(path) %>% rename(run = set) %>%
  mutate(total=in_target_TRUE + in_target_FALSE, 
         set = case_when(
           run == 1 ~ "rar",
           run == 2 ~ "filter_non_target_5sp",
           run == 3 ~ "filter_non_target_humm2anc",
           run == 4 ~ "filter_non_target_humm2anc_honeyanc"
           )) %>%
  select(version, set, gene, count=in_target_TRUE, total)
unique(orig$set)

#read perms
orig_perm<-readRDS(paste0(path_to_data, "/original_galgal6.robj"))
unique(orig_perm$set)

#compute P-value as 1-ecdf(real-1): ecdf(real-1) is prob (X <= (x-1)), e.g. prob(X < x), 1-that is prob (X >= x)

orig_merge <- full_join(orig, orig_perm, 
                        by=c("version" = "version", "set" = "set", "gene" = "gene")) %>% rowwise %>% 
  mutate(pval = 1-ecdf_gene(count-1)) %>% ungroup %>% group_by(version, set) %>% 
  mutate(pval = ifelse(pval == 0, 1e-04, pval)) 

#plot
exclude_nonsig <- function(DF) {
  DF %>% filter(sig_class != "none") %>%
    filter(rar_count > 5 ) # maggie
}

spread_n <- function(df, key, value) {
  # quote key
  keyq <- rlang::enquo(key)
  # break value vector into quotes
  valueq <- rlang::enquo(value)
  s <- rlang::quos(!!valueq)
  df %>% gather(variable, value, !!!s) %>%
    unite(temp, !!keyq, variable) %>%
    spread(temp, value)
}

#clean up merged dataset for plotting
orig_merge_plot <- orig_merge %>%
  group_by(set) %>% 
  filter(version == "basic", count > 0) %>%
  mutate(qval = p.adjust(pval, "fdr")) %>%
  separate(gene, into=c("useless", "sym"), sep="-") %>%
  select(sym, set, cnee_total = total, count, qval) %>%
  ungroup %>%
  spread_n(set, c(qval,count))

path = paste0(project_path, "NumCnee_byGene_", Sys.Date(), ".tsv")       # all data
# write_tsv(orig_merge_plot, path)

#change color scheme, change plotting character
orig_merge_plot_m = orig_merge_plot %>% 
  mutate(sig_class = case_when(
    filter_non_target_5sp_qval <= 0.05 ~ "conv_5sp",
    filter_non_target_humm2anc_honeyanc_qval <= 0.05 ~ "conv_humm2anc_honeyanc",
    filter_non_target_humm2anc_qval <= 0.05 ~ "conv_humm2anc",
    rar_qval <= 0.05 ~ "accelerated",
    TRUE ~ "none"
  )) 

# orig_merge_plot_m = orig_merge_plot %>% 
#   mutate(sig_class = case_when(
#     filter_non_target_humm2anc_qval <= 0.05 ~ "conv_humm2anc",
#     filter_non_target_5sp_qval <= 0.05 ~ "conv_5sp",
#     rar_qval <= 0.05 ~ "accelerated",
#     TRUE ~ "none"
#   )) 
# unique(orig_merge_plot_m$sig_class)

orig_merge_plot_m$sig_class = factor(orig_merge_plot_m$sig_class, 
                                     levels = c("conv_5sp", "conv_humm2anc", "conv_humm2anc_honeyanc",
                                                "accelerated", "none"))
colpal = brewer.pal(4, 'Dark2')
col = c(colpal[4], colpal[2], colpal[3], colpal[1], 'grey70' )

p = orig_merge_plot_m %>%
  ggplot(aes(x=rar_count, y=cnee_total, label=sym)) + 
  geom_jitter(aes(col=sig_class), shape=16, alpha = 0.8) +
  scale_y_log10() + 
  scale_color_manual(values = col) +
  labs(x="Number of convergent accelerated CNEEs near gene", 
       y="Total number of CNEEs near gene", color="Signficantly enriched?") +
  geom_text(data=exclude_nonsig, nudge_x = .2, show.legend=FALSE, 
            size = FontSize, hjust = 0, color = 'grey50') + 
  scale_x_continuous(limits = c(NA, 30), breaks=seq(0, 30, by = 5)) + 
  theme_m + theme(legend.position = c(0.8, 0.2)) ; p

graph_path = paste0(project_path, "NumCnee_byGene_", Sys.Date(), ".svg")       # all data
svg(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*2, pointsize = AxisTxFontSizeSize, onefile = TRUE)
# png(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*0.625, unit = "in", res = 600)
p
dev.off()
openfile(graph_path)


# input
# filter_non_target_5sp: 1
# filter_non_target_humm2anc: 385
# filter_non_target_humm2: 540
# filter_non_target_humm2anc_triMol: 32
# filter_non_target_humm2anc_honeyanc: 16
# filter_non_target_honey2anc_triMol: 2
