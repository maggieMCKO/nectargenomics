#### new perms using shuf #### 
library(RColorBrewer)
library(ggrepel)
library(tidyverse)

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
                panel.border = element_blank(), # element_rect(color = "black", fill= NA, size =.2), 
                axis.line = element_line(color = "black", size = 0.5),
                axis.title.x = element_text(colour = "black", size = AxisTitleFontSizeSize, hjust = 0.5),
                axis.title.y = element_text(colour = "black", size = AxisTitleFontSizeSize),
                axis.text.x = element_text(colour = "black", size = AxisTxFontSizeSize),
                axis.text.y = element_text(colour = "black", size = AxisTxFontSizeSize, hjust = 0),
                axis.ticks = element_line(colour = "black", size = 0.05),
                strip.text = element_text(colour = "black", size = AxisTxFontSizeSize),
                legend.background = element_blank(),
                legend.key = element_blank(),
                legend.key.size = unit(8, "pt"),
                legend.title = element_blank(),
                legend.text = element_text(colour = "black", size = AxisTxFontSizeSize),
                legend.position = "right")

project_path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis_2021/05_perm_phyloacc/1_prep_new/postPhyloAcc/gene/")
target_files = list.files(project_path, ".bed", full.names = T)

out = lapply(target_files, function(s){
  # s = target_files[1]
  print(paste0("current: ", s))
  tmp_name = basename(s)
  data <- read_delim(s, delim = "\t", col_names = F)
  
  data.clean = data %>%
    mutate(X4 = gsub("gene-", "", X4)) 
  
  # calculating pVal this way is off by an order of magnitude .. sum(X7:X1006) isn't working the way I want it to, but count() and n() don't work either and sumRows() will add the actual numbers (e.g. if perm = 2 and is gt obs, then will be summed as 2 perms gt obs instead of one perm)
  pVal <- data.clean %>%  
    rowwise() %>%
    mutate(pVal = sum((X8:X1007 >= X6)+1)/(1000+1)) %>%
    dplyr::select(X1:X7, pVal) %>%
    rename(chr = X1, start= X2, end = X3, gene = X4, biotype = X5, accel = X6, total = X7)
  
  # try calculating pVal this way instead - set logical TRUE for count > obs, convert logical to numerical, then sum per gene 
  t <- data.clean %>% 
    dplyr::select(X4, X5, X6, X7, X8:X1007) %>% 
    dplyr::rename(gene = X4, biotype = X5, obs = X6, total = X7) %>% 
    pivot_longer(cols = starts_with("X"), names_to = "perms", values_to = "count") %>%
    mutate(gt = count >= obs)
  cols <- sapply(t, is.logical) 
  t[,cols] <- lapply(t[,cols], as.numeric)
  tpval <- t %>%
    group_by(gene) %>%
    mutate(sum = sum(gt) + 1,
           pVal = sum/1001)
  tidyt <- tpval %>%
    pivot_wider(names_from = perms, values_from = count) %>%
    dplyr::select(gene, biotype, obs, total, pVal) %>%
    filter(obs > 0) %>%
    distinct() 
  
  # Adjust p old school
  pv <- data.frame(adjP = p.adjust(tidyt$pVal, method = "fdr"), pVal = tidyt$pVal)
  # plot to make sure adjustment worked (ie/ not a 1:1 line)
  plot(-log10(pv$adjP), -log10(pv$pVal))
  # only keep the adjusted p-values
  pv <- pv %>% dplyr::select(-c(pVal))
  adjP <- bind_cols(tidyt, pv) %>%
    mutate(sig_class = case_when(
      adjP <= 0.05 ~ "Enriched for accelerated regions (5% FDR)",
      adjP > 0.05 ~ "Not enriched"))
  
  adjP$sig_class = factor(adjP$sig_class, levels = c("Not enriched", "Enriched for accelerated regions (5% FDR)"))

  
  p = ggplot(adjP, aes(x = obs, y = total, col = sig_class, label = gene)) +
    theme_classic() +
    scale_y_log10() +
    geom_jitter(shape = 16, width = 0.2, show.legend = F, size = 2) +
    scale_colour_brewer(palette = "Dark2") +
    scale_alpha_discrete(range = c(0.5, 0.9)) +
    labs(x = "Number of accelerated CNEEs near gene", y = "Total number of CNEEs near gene", color = "Significance") +
    # scale_x_continuous(breaks = c(1,2,3,4,5)) +
    guides( alpha = "none") +
    geom_text_repel(data = subset(adjP, adjP <= 0.05),
                     aes(label = gene),
                     box.padding = 0.5,
                     point.padding = 0.5,
                     segment.size = 0.1,
                     force = 100,
                     segment.colour = 'grey50') +
    theme(legend.position = "bottom"); p
  
  tmp_name_mod = gsub("\\.bed", paste0("_",  Sys.Date(),  ".pdf"), tmp_name)
  graph_path = paste0(project_path,  tmp_name_mod); graph_path
  pdf(graph_path, width = 6.7, height = 4.2, 
      pointsize = AxisTxFontSizeSize)
  print(p)
  dev.off()
  return(adjP)
})


Ncol = length(target_files)
out_l = lapply(1:Ncol, function(i){
  out[[i]] %>% bind_cols( basename(target_files[i]) )
}) %>% bind_rows() %>%
  mutate(ori = gsub("\\.counts\\.bed", "", ...8),
         ori = gsub("cnee_perms_", "", ori)) %>%
  dplyr::select(-...8)

path = paste0(project_path,  "FullGeneList.tsv")
# write_tsv(out_l, path, na = "")

geneList = out_l %>%
  filter(adjP < 0.05) 
path = paste0(project_path,  "SigGeneList.tsv")
# write_tsv(geneList, path, na = "")

geneList_prot = geneList %>%
  filter(biotype == "protein_coding")
path = paste0(project_path,  "SigGeneList_proteincoding.tsv")
# write_tsv(geneList_prot, path, na = "")

## bin
geneList_sum = geneList %>% 
  group_by(ori) %>%
  tally() %>%
  mutate(ori = fct_reorder(ori, desc(n)))

p_bins = ggplot(geneList_sum, aes(x = ori, y = n)) +
  geom_col(fill = 'steelblue') +
  geom_text(aes(label = n), vjust = -1, size = FontSize) +
  scale_x_discrete(expand = c(0,0,0.05,0)) +
  scale_y_continuous(expand = c(0,0,0.05,0), "Count") +
  theme_m + 
  theme(axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()); p_bins

graph_path = paste0(project_path,  "gene_count.pdf"); graph_path
pdf(graph_path, width = 4.2, height = 3.35, pointsize = AxisTxFontSizeSize)
p_bins
dev.off()

#### upset =====

library(UpSetR)

present = function(s){ifelse(!is.na(s), 1, 0)}

geneList_wide = geneList %>%
  ## to wide
  pivot_wider(id_cols = gene, names_from = ori, values_from = c(adjP), values_fn = present, values_fill = 0) %>%
  as.data.frame()

p = upset(geneList_wide,
          # sets = c("swifts", "falcons", "lyrebirds", "passerides"), keep.order = TRUE,
          nsets = 3,
          main.bar.color = "#56B4E9",
          matrix.color = "#56B4E9",
          sets.bar.color = "#56B4E9",
          # order.by = c("freq", "degree"),
          empty.intersections = "on"); p


graph_path = paste0(project_path, "Gene_upset_", Sys.Date(),".pdf"); graph_path
pdf(graph_path, width = 6.7, height = 3.35, pointsize = AxisTxFontSizeSize)
p
dev.off()
