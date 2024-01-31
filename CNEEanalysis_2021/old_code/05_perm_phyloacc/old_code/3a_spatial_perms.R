library(tidyverse)
library(qqman)

# ref: https://github.com/sjswuitchik/duck_comp_gen/blob/master/03_cnee_analyses/05a_cnees_window.R

project_path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis_2021/05_perm_phyloacc/1_prep_new/postPhyloAcc/spatial/")
target_files = list.files(project_path, "acc", full.names = T); target_files

df = tibble(sets = c("hummingbirds", "parrots", "honeyeaters_pardalote", "sunbirds_flowerpecker"),
            n_perms = c(1691, 1125, 655, 481))

# func ####
out = lapply(target_files, function(s){
  # s = target_files[1]
  print(paste0("current: ", s))
  tmp_name = basename(s)
  tmp_name = gsub("window.acc.cnees_", "", tmp_name)
  tmp_name = gsub(".bed", "", tmp_name)
  
  # load in accelerated CNEEs by window
  dataa <- read_delim(s, delim = "\t", col_names = c("chr", "start", "end", "acc_ori")) # 211978
  
  acc.clean.window = dataa %>%
    # convert remaining CNEE names to a binary measure of occurrence (1)
    mutate(acc = ifelse(grepl("CNEE", acc_ori), 1, 0)) %>%
    select(-acc_ori) %>%
    # summarise by window
    group_by(chr, start, end) %>%
    summarise(acc.cnee = sum(acc)) # 211848
  sum(acc.clean.window$acc.cnee) # 821
  
  # load in total list of CNEEs by window
  path = paste0(project_path, "../window.cnees.bed")
  data_all = read_delim(path, delim = "\t", col_names = c("chr", "start", "end", "cnee_ori")) # 638942

  all.clean.window = data_all %>%
    # convert CNEE names to a binary measure of occurrence (1)
    mutate(cnee = ifelse(grepl("CNEE", cnee_ori), 1, 0)) %>%
    select(-cnee_ori) %>%
    # summarise by window
    group_by(chr, start, end) %>%
    summarise(total.cnee = sum(cnee)) # 211848
  sum(all.clean.window$total.cnee) # 439254
    
  # combine, clean, and calculate proportion of accel'd CNEEs per window
  data.window <- full_join(all.clean.window, acc.clean.window, by = c("chr", "start", "end")) %>%
    filter(total.cnee > 0) # 1: 12160
  final.data.window <- data.window %>%
    mutate(prop = acc.cnee/total.cnee) 
  
  # binomial test function
  n_tmp = df$n_perms[match(tmp_name, df$sets)]
  t = 363747
  
  bt <- function(x, n, p = n_tmp/t) {
    binom.test(x, n, n_tmp/t, alternative = "greater", conf.level = 0.95)$p.value
  }
  
  # add binomial p-values to table
  final.data.window$pVal <- mapply(bt, final.data.window$acc.cnee, final.data.window$total.cnee)
  
  # adjust p-values for multiple comparisons
  # don't use mutate for this, do it old school
  pv <- data.frame(adjustP = p.adjust(final.data.window$pVal, method = "fdr"), pVal = final.data.window$pVal)
  # plot to make sure adjustment worked (ie/ not a 1:1 line)
  plot(-log10(pv$adjustP), -log10(pv$pVal))
  # only keep the adjusted p-values
  pv <- pv %>% select(-c(pVal))
  
  adjP <- bind_cols(final.data.window, pv)
  
  
  test <- left_join(dataa, adjP, by = c("chr" = "chr", "start" = "start", "end" = "end")) %>%
    na.omit()  %>% 
    mutate(chr_m = gsub("chr", "", chr)) %>% 
    mutate(chr_m = gsub("9_.*", 9.2, chr_m)) %>% 
    mutate(chr_m = gsub("W", 34, chr_m)) %>%
    mutate(chr_m = gsub("Z", 35, chr_m)) %>%
    mutate(chr_m = as.numeric(chr_m))# %>%
    # mutate(lab = ifelse(adjustP < 0.05, acc_ori, NA)) %>%
    # left_join(anno, by = c("acc_ori" = "ID")) 
  path = paste0(project_path, tmp_name, "_testVal.tsv")
  # write_tsv(test, path, col_names = T, na= "")
  # unique(test$chr_m)
  
  test_sub = test %>% select(-acc_ori) %>% distinct()
  
  qq(test$adjustP)
  
  # pval distribution
  p = ggplot(adjP) + 
    geom_density(aes(x= pVal)) +
    ggtitle(tmp_name) +
    theme_classic() + theme(plot.title = element_text(size = AxisTitleFontSizeSize, hjust = 0.5) ); p
  graph_path = paste0(project_path, tmp_name, "_pVal_dist_", Sys.Date(), ".pdf")
  pdf(graph_path, width = 4.2, height = 3.35, pointsize = AxisTxFontSizeSize, onefile = TRUE)
  print(p)
  dev.off()
  
  
  p = ggplot(adjP) + 
    geom_density(aes(x= adjustP), color = "red") +
    # scale_x_continuous(limits = c(NA, 0.1)) +
    # scale_y_continuous(limits = c(NA, 10)) +
    ggtitle(tmp_name) +
    theme_classic() + theme(plot.title = element_text(size = AxisTitleFontSizeSize, hjust = 0.5) ); p
  graph_path = paste0(project_path, tmp_name, "_padj_dist_", Sys.Date(), ".pdf")
  pdf(graph_path, width = 4.2, height = 3.35, pointsize = AxisTxFontSizeSize, onefile = TRUE)
  print(p)
  dev.off()
  
  
  # # FDR 1%
  # gwide <- test %>% filter( adjustP < 0.01) %>%
  ##   slice_max(order_by = adjustP, n = 1, with_ties = F ) %>%
  #  summarize(maxpval= max(pVal)) %>% pull(maxpval)
  # 
  # # FDR 5%
  # ss <- test %>% filter( adjustP < 0.05) %>%
  ##   slice_max(order_by = adjustP, n = 1, with_ties = F ) %>%
  #   summarize(maxpval= max(pVal)) %>% pull(maxpval)
  
  adjP_sub = adjP %>% filter(adjustP < 0.05) %>%
    left_join(dataa, by = c("chr" = "chr", "start" = "start", "end" = "end")) %>%
    left_join(anno, by = c("acc_ori" = "ID"))
  
  graph_path = paste0(project_path, tmp_name, "_manhattan_", "pval", "_", Sys.Date(), ".pdf")
  pdf(graph_path, width = 7.7, height = 3.35, pointsize = AxisTxFontSizeSize, onefile = TRUE)
  manhattan(test, chr="chr_m", 
            bp = "start", snp = "acc_ori", p = "pVal", col = c("grey", "skyblue"), 
            # suggestiveline=FALSE, genomewideline = FALSE
            # annotatePval = 1e-05,
            # annotateTop = TRUE
            # chrlabs = "chr", 
            # highlight = test$lab %>% na.omit(),
            # suggestiveline=-log10(ss), genomewideline = -log10(gwide), 
            # ylim = c(0,8)
  ) # note: W = 34, Z = 35
  dev.off()
  
  
  if(nrow(adjP_sub) > 0){
    
    path = paste0(project_path, tmp_name, "_final.cnees.window.sig.tsv")
    write_tsv(adjP_sub, path, col_names = T, na = "")
    
    path = paste0(project_path, tmp_name, "_final.cnees.window.sig.bed")
    adjP_sub2 = adjP_sub %>% select(chr, start, end) %>% distinct()
    write_tsv(adjP_sub2, path, col_names = F, na = "")
    
  }

})

###  run 3b_anno_peaks.sh to anno windows

# combind results ####
target_files2 = list.files(project_path, "sig_intersect.bed", full.names = T); target_files2

# func ####
out2 = lapply(target_files2, function(s){
  # s = target_files2[1]
  print(paste0("current: ", s))
  tmp_name = basename(s)
  tmp_name = gsub("_final.cnees.window.sig_intersect.bed", "", tmp_name)
  
  # load in accelerated CNEEs by window
  dataa <- read_tsv(s, col_names = c("chr", "start", "end", "chr_gene", "start_gene", "end_gene", "gene", "biotype")) %>%
    bind_cols('group' = tmp_name) %>% 
    filter(biotype == "protein_coding") %>%
    select(group, chr, start, end, gene, biotype)
  
}) %>% bind_rows()

path = paste0(project_path, "AccCneeEnrichedWindows_proteincoding_padj0.05.tsv")
write_tsv(out2, path, col_names = T, na = "")
