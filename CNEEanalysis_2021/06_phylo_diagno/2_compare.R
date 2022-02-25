library("tidyverse")
library("ggtree")
library("grid")
library("gridExtra")
library("RColorBrewer")
# library("gplots")
# library("gtools")

openfile = function(filepath){
  sysinf = Sys.info()
  if (!is.null(sysinf)){
    os = sysinf['sysname']
    if (os == 'Darwin')
      os = "osx"
  } else { ## mystery machine
    os = .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os = "osx"
    if (grepl("linux-gnu", R.version$os))
      os = "linux"
  }
  type = tolower(os);
  switch(os, Linux = system(paste0("xdg-open ", filepath)),
         Windows = system(paste0("open \"", filepath, "\"")),
         osx = system(paste0("open \"", filepath, "\""))
  )}

FontSize = 1.5
AxisTxFontSizeSize = 5
AxisTitleFontSizeSize = 6
Factor_mmtoin = 0.0393701
Width_HalfCol = 85*Factor_mmtoin 
theme_m = theme(panel.background = element_blank(),
                plot.background = element_blank(), # default is white
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
                legend.title = element_blank(),
                legend.text = element_text(colour = "black", size = AxisTxFontSizeSize),
                legend.position = "bottom")

# load files ====
project_path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis_2021/06_phylo_diagno/")

# run1 output
path = paste0(project_path, "../../CNEEanalysis_2021_older/03_phyloacc/1_run_phyloacc/rate_postZ_M2_combined_Target_NeFr40.txt")
run1_rates = read_tsv(path) %>% filter(!is.na(n_rate)) # 310886

path = paste0(project_path, "../../CNEEanalysis_2021_older/03_phyloacc/1_run_phyloacc/phyloacc_score_postZ_filtered_Target_NeFr40.tsv.gz")
run1_df = read_tsv(path) # 310886

run1_df2 = run1_rates[, 1:6] %>% full_join(run1_df)

# run2 output
path = paste0(project_path, "/../03_phyloacc/oldway/output/rate_postZ_M2_combined_Target_NeFr40.txt")
run2_rates = read_tsv(path) %>% filter(!is.na(n_rate)) # 346079

path = paste0(project_path, "/../03_phyloacc/oldway/4_post_phyloacc/phyloacc_score_postZ_filtered_Target_NeFr40.tsv.gz")
run2_df = read_tsv(path) # 346079
# 363658 rows in elem_lik
# conserved set: 363747

run2_df2 = run2_rates[, 1:6] %>% full_join(run2_df)

# phyloP
path = paste0(project_path, "/../01b_cnee_prep/phylop_allChr_lrt_con.tsv")
# path = paste0(project_path, "/../01b_cnee_prep/phylop_allChr_lrt_conacc.tsv")
phylop = read_tsv(path) %>%
  rename("chr" = "#chr") %>%
  filter(!grepl("#chr", chr)) %>%
  rename("id_run2" = "name", "scale_run2" = "scale", 
         "lnlratio_run2" = "lnlratio", "pval_run2" = "pval") %>%
  mutate(padj_run2 = p.adjust(pval_run2, "fdr"), 
         pval_fix_run2 = ifelse(pval_run2 == 0, 10^-10, pval_run2),
         log10pval_run2 = -log10(pval_fix_run2)) # 378215
range(phylop$scale_run2,na.rm = T)

# conversion matrix
path = paste0(project_path, "comparing_cleaned_cnees.bed")
transMatrix = read_tsv(path, col_names = c("chr_run1", "start_run1", "end_run1", "id_run1",
                                           "chr_run2", "start_run2", "end_run2", "id_run2")) 
transMatrix2 = transMatrix %>%
  mutate(len_run1 = end_run1 - start_run1,
         len_run2 = end_run2 - start_run2) %>%
  left_join(run1_df2[, 2:15], by = c("id_run1" = "ID"))%>%
  left_join(run2_df2[, 2:15], by = c("id_run2" = "ID"), suffix = c("_run1", "_run2")) %>%
  mutate(logBF3_run1 = loglik_Full_run1 - loglik_Null_run1,
         logBF3_run2 = loglik_Full_run2 - loglik_Null_run2) %>%
  left_join(phylop)
range(transMatrix2$scale_run2,na.rm = T)
length(setdiff(phylop$id_run2, transMatrix2$id_run2)) # 22045 


# # interface
# # final cnees
# path = paste0("/Users/beautibabi/ownCloud/Baldwin/Projects/Nectar/Seq_Data/globus/CNEEanalysis_2021/01b_cnee_prep/bed_ouputs/galGal6_final_conserved_CNEEs.bed")
# final_bed = read_tsv(path, col_names = c('chr', 'start', 'end', 'id'))
# 
# path = paste0(project_path, "../03_phyloacc/phyloacc-out-02-14-2022.09-46-46/results/id-key.txt")
# run2i_id = read_tsv(path, col_names = F) %>%
#   rename("Locus ID" = "X1", "ind" = "X2") # 361779
# run2i_id$id = final_bed$id[run2i_id$ind]
# 
# path = paste0(project_path, "../03_phyloacc/phyloacc-out-02-14-2022.09-46-46/results/elem_lik.txt")
# run2i_lik = read_tsv(path) %>% 
#   left_join(run2i_id %>% select(`Locus ID`, id)) # 361779
# 
# path = paste0(project_path, "../03_phyloacc/phyloacc-out-02-14-2022.09-46-46/results/rate_postZ_M2.txt") # pp
# run2i_pp = read_tsv(path)  
# run2i_pp2 = run2i_pp %>% select(1:6, matches("_3")) %>% 
#   left_join(run2i_id %>% select(`Locus ID`, id)) # 350578
# 
# transMatrixI = run2i_pp2 %>% select(1:6, id) %>%
#   full_join(run2i_lik %>% select(-`Locus ID`)) %>%
#   full_join(run2_df2[, 2:15], by = c("id" = "ID"), suffix = c("_run2i", "_run2")) %>%
#   mutate(logBF3_run2i = loglik_Full_run2i - loglik_Null_run2i,
#        logBF3_run2 = loglik_Full_run2 - loglik_Null_run2)

# # plot0: run2: interface vs old ====
# # n_rate
# corv = signif(cor(transMatrixI$n_rate_run2i, transMatrixI$n_rate_run2, "complete.obs"), 2); corv # 0.56
# p_n_rate_cori = ggplot(transMatrixI, aes(x = n_rate_run2i, y = n_rate_run2)) +
#   geom_point(size = 0.2, alpha = 0.5) + 
#   geom_abline(intercept = 0, slope = 1, color = 'grey50', size = 0.2, linetype = 2) +
#   geom_abline(intercept = 0, slope = corv, color = 'red') +
#   annotate('text', x = 10, y = 16, 
#            label = paste0("cor: ", corv),
#            size = FontSize, hjust = 0, color = 'red') +
#   
#   theme_m; p_n_rate_cori


# 1. plot1: run1 vs run2 ====
# n_rate
corv = signif(cor(transMatrix2$n_rate_run1, transMatrix2$n_rate_run2, "complete.obs"), 2); corv # 0.56
p_n_rate_cor = ggplot(transMatrix2, aes(x = n_rate_run1, y = n_rate_run2)) +
  geom_point(size = 0.2, alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, color = 'grey50', size = 0.2, linetype = 2) +
  geom_abline(intercept = 0, slope = corv, color = 'red') +
  annotate('text', x = 10, y = 16, 
           label = paste0("cor: ", corv),
           size = FontSize, hjust = 0, color = 'red') +
  
  theme_m; p_n_rate_cor

# c_rate
corv = signif(cor(transMatrix2$c_rate_run1, transMatrix2$c_rate_run2, "complete.obs"), 2); corv # 0.83
p_c_rate_cor = ggplot(transMatrix2, aes(x = c_rate_run1, y = c_rate_run2)) +
  geom_point(size = 0.2, alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, color = 'grey50', size = 0.2, linetype = 2) +
  geom_abline(intercept = 0, slope = corv, color = 'red') +
  annotate('text', x = 0.65, y = 0.45, 
           label = paste0("cor: ", corv),
           size = FontSize, hjust = 1, color = 'red') +
  theme_m; p_c_rate_cor

# len cor
corv = signif(cor(transMatrix2$len_run1, transMatrix2$len_run2), 2); corv # 0.81
p_len_cor = ggplot(transMatrix2, aes(x = len_run1, y = len_run2)) +
  geom_point(size = 0.2, alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, color = 'grey50', size = 0.2, linetype = 2) +
  geom_abline(intercept = 0, slope = corv, color = 'red') +
  annotate('text', x = 5000, y = 3000, 
           label = paste0("cor: ", corv),
           size = FontSize, hjust = 0, color = 'red') +
  theme_m; p_len_cor

p_len_cor_zoomin = p_len_cor + coord_cartesian(xlim = c(NA, 5000)); p_len_cor_zoomin

# loglik_Full
corv = signif(cor(transMatrix2$loglik_Full_run1, transMatrix2$loglik_Full_run2, "complete.obs"), 2); corv # 0.93
p_llM2_cor = ggplot(transMatrix2, aes(x = loglik_Full_run1, y = loglik_Full_run2)) +
  geom_point(size = 0.2, alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, color = 'grey50', size = 0.2, linetype = 2) +
  geom_abline(intercept = 0, slope = corv, color = 'red') +
  annotate('text', x = -1500, y = -8000, 
           label = paste0("cor: ", corv),
           size = FontSize, hjust = 0, color = 'red') +
  theme_m; p_llM2_cor

# loglik_Null
corv = signif(cor(transMatrix2$loglik_Null_run1, transMatrix2$loglik_Null_run2, "complete.obs"), 2); corv # 0.93
p_llM0_cor = ggplot(transMatrix2, aes(x = loglik_Null_run1, y = loglik_Null_run2)) +
  geom_point(size = 0.2, alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, color = 'grey50', size = 0.2, linetype = 2) +
  geom_abline(intercept = 0, slope = corv, color = 'red') +
  annotate('text', x = -2000, y = -8000, 
           label = paste0("cor: ", corv),
           size = FontSize, hjust = 0, color = 'red') +
  theme_m; p_llM0_cor

# loglik_Acc
corv = signif(cor(transMatrix2$loglik_Acc_run1, transMatrix2$loglik_Acc_run2, "complete.obs"), 2); corv # 0.93
p_llM1_cor = ggplot(transMatrix2, aes(x = loglik_Acc_run1, y = loglik_Acc_run2)) +
  geom_point(size = 0.2, alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, color = 'grey50', size = 0.2, linetype = 2) +
  geom_abline(intercept = 0, slope = corv, color = 'red') +
  annotate('text', x = -2000, y = -8000, 
           label = paste0("cor: ", corv),
           size = FontSize, hjust = 0, color = 'red') +
  theme_m; p_llM1_cor

# logBF1
corv = signif(cor(transMatrix2$logBF1_run1, transMatrix2$logBF1_run2, "complete.obs"), 2); corv # 0.064
p_bf1_cor = ggplot(transMatrix2, aes(x = logBF1_run1, y = logBF1_run2)) +
  geom_point(size = 0.2, alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, color = 'grey50', size = 0.2, linetype = 2) +
  geom_abline(intercept = 0, slope = corv, color = 'red') +
  annotate('text', x = 200, y = 150, 
           label = paste0("cor: ", corv),
           size = FontSize, hjust = 0, color = 'red') +
  theme_m; p_bf1_cor

# logBF2
corv = signif(cor(transMatrix2$logBF2_run1, transMatrix2$logBF2_run2, "complete.obs"), 2); corv # 0.56
p_bf2_cor = ggplot(transMatrix2, aes(x = logBF2_run1, y = logBF2_run2)) +
  geom_point(size = 0.2, alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, color = 'grey50', size = 0.2, linetype = 2) +
  geom_abline(intercept = 0, slope = corv, color = 'red') +
  annotate('text', x = -250, y = -80, 
           label = paste0("cor: ", corv),
           size = FontSize, hjust = 0, color = 'red') +
  theme_m; p_bf2_cor

# logBF3
corv = signif(cor(transMatrix2$logBF3_run1, transMatrix2$logBF3_run2, "complete.obs"), 2); corv # 0.47
p_bf3_cor = ggplot(transMatrix2, aes(x = logBF3_run1, y = logBF3_run2)) +
  geom_point(size = 0.2, alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, color = 'grey50', size = 0.2, linetype = 2) +
  geom_abline(intercept = 0, slope = corv, color = 'red') +
  annotate('text', x = 400, y = 300, 
           label = paste0("cor: ", corv),
           size = FontSize, hjust = 0, color = 'red') +
  theme_m; p_bf3_cor

# combine
p = arrangeGrob(grobs = list(p_len_cor, p_c_rate_cor, p_n_rate_cor,
                             p_bf1_cor, p_bf2_cor, p_bf3_cor, 
                             p_llM0_cor, p_llM1_cor, p_llM2_cor), ncol = 3, )

graph_path = paste0(project_path, "comparingRuns_",  Sys.Date(),  ".png"); graph_path
png(graph_path, width = Width_HalfCol*1.5, height = Width_HalfCol*1.5, units = 'in', res = 300,
    pointsize = AxisTxFontSizeSize)
grid.draw(p)
dev.off()
openfile(file = graph_path)

# 2. plot2: run2: phyloacc vs phylop ====
# run2: logBF3 vs phylop pval
# corv = signif(cor(transMatrix2$logBF3_run1, transMatrix2$logBF3_run2, "complete.obs"), 2); corv # 1

# p_bf3_pval_cor = ggplot(transMatrix2, aes(x = log10pval_run2, y = logBF3_run2)) +
#   geom_point(size = 0.2, alpha = 0.5) + 
#   geom_abline(intercept = 0, slope = 1, color = 'red') +
#   theme_m; p_bf3_pval_cor
# 
# p_bf3_pval_cor_zoomin = p_bf3_pval_cor + coord_cartesian(xlim = c(NA, 6)); p_bf3_pval_cor_zoomin


p_bf3_pvalraw_cor = ggplot(transMatrix2, aes(x = -log10(pval_run2), y = logBF3_run2)) +
  geom_point(size = 0.2, alpha = 0.5) + 
  # geom_abline(intercept = 0, slope = 1, color = 'red') +
  theme_m; p_bf3_pvalraw_cor

# p_bf3_pvalraw_cor_zoomin = p_bf3_pvalraw_cor + coord_cartesian(ylim = c(NA, 0)); p_bf3_pvalraw_cor_zoomin

# p_bf3_pvalraw_cor = ggplot(transMatrix2, aes(x = pval_run2, y = logBF3_run2)) +
#   geom_point(size = 0.2, alpha = 0.5) + 
#   # geom_abline(intercept = 0, slope = 1, color = 'red') +
#   theme_m; p_bf3_pvalraw_cor
# 
# p_bf3_pvalraw_cor2 = ggplot(transMatrix2, aes(x = log10(pval_run2), y = logBF3_run2)) +
#   geom_point(size = 0.2, alpha = 0.5) + 
#   # geom_abline(intercept = 0, slope = 1, color = 'red') +
#   theme_m; p_bf3_pvalraw_cor2

# run2: c_rate vs phylop scale
corv = signif(cor(transMatrix2$scale_run2, transMatrix2$c_rate_run2, "complete.obs"), 2); corv # 0.56
p_c_rate_cor = ggplot(transMatrix2, aes(x = scale_run2, y = c_rate_run2)) +
  geom_point(size = 0.2, alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, color = 'grey50', size = 0.2, linetype = 2) +
  geom_abline(intercept = 0, slope = corv, color = 'red') +
  annotate('text', x = 0.7, y = 0.6, 
           label = paste0("cor: ", corv),
           size = FontSize, hjust = 0, color = 'red') +
  theme_m; p_c_rate_cor

# run2: n_rate vs phylop scale
corv = signif(cor(transMatrix2$scale_run2, transMatrix2$n_rate_run2, "complete.obs"), 2); corv # 0.049
p_n_rate_cor = ggplot(transMatrix2, aes(x = scale_run2, y = n_rate_run2)) +
  geom_point(size = 0.2, alpha = 0.5) + 
  # geom_abline(intercept = 0, slope = 1, color = 'grey50', size = 0.2, linetype = 2) +
  # geom_abline(intercept = 0, slope = corv, color = 'red') +
  # annotate('text', x = 0.7, y = 0.6, 
  #          label = paste0("cor: ", corv),
  #          size = FontSize, hjust = 0, color = 'red') +
  theme_m; p_n_rate_cor

p = arrangeGrob(grobs = list(p_bf3_pvalraw_cor, p_c_rate_cor), nrow = 1)

# graph_path = paste0(project_path, "comparing_run2_phyloacc_phyloP_",  Sys.Date(),  ".png"); graph_path
graph_path = paste0(project_path, "comparing_run2_phyloacc_phyloP_conacc_",  Sys.Date(),  ".png"); graph_path
png(graph_path, width = Width_HalfCol*1.2, height = Width_HalfCol*0.6, units = 'in', res = 300,
    pointsize = AxisTxFontSizeSize)
grid.draw(p)
dev.off()
openfile(file = graph_path)

# 3. tree ====
x = "((((((((((((((((((HLmelMel1:0.0521152,HLserCan1:0.0340235)HLmelMel1-HLserCan1:0.020135,HLtaeGut4:0.0443394)HLmelMel1-HLtaeGut4:0.0119299,((HLlepAsp1:0.0226641,HLcinPul1:0.022812)HLlepAsp1-HLcinPul1:0.0106772,HLdicExi1:0.0421218)HLlepAsp1-HLdicExi1:0.0108789)HLmelMel1-HLlepAsp1:0.00596355,ficAlb2:0.0675878)HLmelMel1-ficAlb2:0.00480446,(HLparMaj1:0.0142069,pseHum1:0.0128063)HLparMaj1-pseHum1:0.0380282)HLmelMel1-HLparMaj1:0.0161352,HLcorCor3:0.0455935)HLmelMel1-HLcorCor3:0.0099441,((HLmalCya1:0.018977,HLmalEle1:0.0156481)HLmalCya1-HLmalEle1:0.0605155,((((HLlicMelCas1:0.0132722,HLlicPen1:0.0133032)HLlicMelCas1-HLlicPen1:0.00760819,HLphyNov1:0.0181313)HLlicMelCas1-HLphyNov1:0.00826567,HLgraPic1:0.0216377)HLlicMelCas1-HLgraPic1:0.033752,(HLacaPus1:0.0409313,HLparPun1:0.0422812)HLacaPus1-HLparPun1:0.00771203)HLlicMelCas1-HLacaPus1:0.00880152)HLmalCya1-HLlicMelCas1:0.0132206)HLmelMel1-HLmalCya1:0.0166602,HLcliRuf1:0.0783201)HLmelMel1-HLcliRuf1:0.0137818,(HLatrCla1:0.0449263,HLmenNov1:0.035249)HLatrCla1-HLmenNov1:0.020087)HLmelMel1-HLatrCla1:0.0318542,(HLfurRuf1:0.0769285,HLempTra1:0.0782382)HLfurRuf1-HLempTra1:0.0284616)HLmelMel1-HLfurRuf1:0.0774837,((((HLamaAes1:0.0255052,HLaraSol1:0.0271174)HLamaAes1-HLaraSol1:0.0140122,(HLtriMol2:0.0312032,HLlorGal1:0.0365652)HLtriMol2-HLlorGal1:0.0102388)HLamaAes1-HLtriMol2:0.00891485,HLnymHol2:0.036825)HLamaAes1-HLnymHol2:0.0162666,HLstrHab1:0.039069)HLamaAes1-HLstrHab1:0.0544873)HLmelMel1-HLamaAes1:0.0148137,((falPer1:0.00269468,falChe1:0.00260786)falPer1-falChe1:0.00628114,HLfalTin1:0.00784217)falPer1-HLfalTin1:0.0739271)HLmelMel1-falPer1:0.00516277,(HLtytAlb2:0.0687773,halLeu1:0.0539185)HLtytAlb2-halLeu1:0.00651593)HLmelMel1-HLtytAlb2:0.00424469,aptFor1:0.0466123)HLmelMel1-aptFor1:0.0058816,(HLcalPug1:0.106203,opiHoa1:0.0790076)HLcalPug1-opiHoa1:0.00896489)HLmelMel1-HLcalPug1:0.00788481,((((HLfloFus1:0.0305673,HLphaSup1:0.0312286)HLfloFus1-HLphaSup1:0.00470972,HLcalAnn5:0.0407605)HLfloFus1-HLcalAnn5:0.0806945,(HLchaPel1:0.0335644,HLapuApu1:0.0236387)HLchaPel1-HLapuApu1:0.0676889)HLfloFus1-HLchaPel1:0.0327126,cucCan1:0.121706)HLfloFus1-cucCan1:0.00632473)HLmelMel1-HLfloFus1:0.0130101,HLcolLiv2:0.113698)HLmelMel1-HLcolLiv2:0.0966404,galGal6:0.0966404)HLmelMel1-galGal6;
"
tree = read.tree(text=x)

x_lim = max(tree$edge.length)*4
FontSize_tree = 3
tree_p = ggtree(tree, layout = "slanted") +
  geom_tiplab( offset = 0, size = FontSize_tree) +
  scale_color_manual(values = c('black', 'blue')) +
  xlim_tree(x_lim) +
  geom_nodelab(size = FontSize, nudge_y = 0.2) +
theme_tree() + theme(legend.position = 'none') ; tree_p
tree_p = tree_p %>% ggtree::rotate(46) ; tree_p


graph_path = paste0(project_path, "tree_",  Sys.Date(),  ".png"); graph_path
png(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*4, units = 'in', res = 300,
    pointsize = AxisTxFontSizeSize)
tree_p
dev.off()
openfile(file = graph_path)

# 4. total number of branches > 0.95 ====
# cal_branch = function(input_df){
#   # input_df = run1_df
#   
#   # keep only posterior probability and CNEE ID
#   input_df = input_df %>% select(!matches("^lo") ) %>% 
#     select(!matches("basic")) %>% select(!matches("filtered")) %>%
#     select(-No.) 
#   
#   # matrix with only posterior probability
#   input_dfm = input_df %>% select(-ID) %>% as.matrix()
#   x = input_dfm >= 0.95
#   sum_row = apply(x, 1, sum)
#   names(sum_row) = input_df$ID
#   sum_row_tb = sum_row %>% as_tibble(rownames = NA) %>% 
#     rownames_to_column(var = 'id') %>% rename('count' = 'value')
# }
# 
# run1_count = cal_branch(run1_df)
# run2_count = cal_branch(run2_df)

cal_only_tips = function(input_df){
  # input_df = run1_df
  
  # keep only posterior probability and CNEE ID
  input_df = input_df %>% select(!matches("^lo") ) %>% 
    select(!matches("basic")) %>% select(!matches("filtered")) %>%
    select(!matches("-")) %>%
    select(-No.) 
  
  # matrix with only posterior probability
  input_dfm = input_df %>% select(-ID) %>% as.matrix()
  x = input_dfm >= 0.95
  sum_row = apply(x, 1, sum)
  names(sum_row) = input_df$ID
  sum_row_tb = sum_row %>% as_tibble(rownames = NA) %>% 
    rownames_to_column(var = 'id') %>% rename('count' = 'value')
}

run1_count = cal_only_tips(run1_df)
run2_count = cal_only_tips(run2_df)

transMatrix3_c = transMatrix %>%
  left_join(run1_count, by = c("id_run1" = "id")) %>%
  left_join(run2_count, by = c("id_run2" = "id"), suffix = c("_run1", "_run2")) %>%
  left_join(phylop) %>%
  left_join(run2_df2[, 2:15], by = c("id_run2" = "ID")) %>%
  mutate(logBF3_run2 = loglik_Full - loglik_Null)

# expectation of acceleration
corv = signif(cor(transMatrix3_c$count_run1, transMatrix3_c$count_run2, "complete.obs"), 2); corv # 0.79
p_count_cor = ggplot(transMatrix3_c, aes(x = count_run1, y = count_run2)) +
  geom_point(size = 0.2, alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, color = 'grey50', size = 0.2, linetype = 2) +
  geom_abline(intercept = 0, slope = corv, color = 'red') +
  annotate('text', x = 10, y = 16, 
           label = paste0("cor: ", corv),
           size = FontSize, hjust = 0, color = 'red') +
  theme_m; p_count_cor

## histogram ====
nbins = 50

# count run1
p_count1 = ggplot(transMatrix3_c, aes(x = count_run1)) +
  geom_histogram(bins = nbins) +
  scale_y_continuous(expand = c(0,0,0.05,0), limits = c(NA, 35000)) +
  # geom_vline(xintercept = 0, color = 'red') +
  # annotate("text", x = 0, y = 20000,
  #          label = "x = 1",
  #          size = FontSize*1.5, hjust = -0.5, color = 'red') +
  ggtitle(paste0("run1")) +
  theme_m; p_count1

# count run2
p_count2 = ggplot(transMatrix3_c, aes(x = count_run2)) +
  geom_histogram(bins = nbins) +
  scale_y_continuous(expand = c(0,0,0.05,0), limits = c(NA, 35000)) +
  # geom_vline(xintercept = 0, color = 'red') +
  # annotate("text", x = 0, y = 20000,
  #          label = "x = 1",
  #          size = FontSize*1.5, hjust = -0.5, color = 'red') +
  ggtitle(paste0("run2")) +
  theme_m; p_count2

p = arrangeGrob(grobs = list(p_count1, p_count2), nrow = 1)

graph_path = paste0(project_path, "comparing_runs_NumBranch_gl95_",  Sys.Date(),  ".png"); graph_path
png(graph_path, width = Width_HalfCol*1, height = Width_HalfCol*0.4, units = 'in', res = 300,
    pointsize = AxisTxFontSizeSize)
grid.draw(p)
dev.off()
openfile(file = graph_path)


## compare with phyloP ====
# count vs pval
p_count2_pvalraw_cor = ggplot(transMatrix3_c, aes(x = pval_run2, y = count_run2)) +
  geom_point(size = 0.2, alpha = 0.5) + 
  # geom_abline(intercept = 0, slope = 1, color = 'red') +
  theme_m; p_count2_pvalraw_cor

# count vs scale
corv = signif(cor(transMatrix3_c$scale_run2, transMatrix3_c$count_run2, "complete.obs"), 2); corv # 0.57
p_count2_scale_cor = ggplot(transMatrix3_c, aes(x = scale_run2, y = count_run2)) +
  geom_point(size = 0.2, alpha = 0.5) + 
  # geom_abline(intercept = 0, slope = 1, color = 'grey50', size = 0.2, linetype = 2) +
  # geom_abline(intercept = 0, slope = corv, color = 'red') +
  annotate('text', x = 0.2, y = 20, 
           label = paste0("cor: ", corv),
           size = FontSize, hjust = 0, color = 'red') +
  theme_m; p_count2_scale_cor

# count vs logBF3
corv = signif(cor(transMatrix3_c$logBF3_run2, transMatrix3_c$count_run2, "complete.obs"), 2); corv # -0.66
p_count2_bf3_cor = ggplot(transMatrix3_c, aes(x = logBF3_run2, y = count_run2)) +
  geom_point(size = 0.2, alpha = 0.5) + 
  # geom_abline(intercept = 0, slope = 1, color = 'grey50', size = 0.2, linetype = 2) +
  # geom_abline(intercept = 0, slope = corv, color = 'red') +
  # annotate('text', x = 100, y = 50, 
  #          label = paste0("cor: ", corv),
  #          size = FontSize, hjust = 0, color = 'red') +
  theme_m ; p_count2_bf3_cor

p = arrangeGrob(grobs = list(p_count2_pvalraw_cor, p_count2_scale_cor, p_count2_bf3_cor), nrow = 1)

graph_path = paste0(project_path, "comparing_run2_phyloacc_phyloP_NumBranch_gl95_",  Sys.Date(),  ".png"); graph_path
# graph_path = paste0(project_path, "comparing_run2_phyloacc_phyloPconacc_NumBranch_gl95_",  Sys.Date(),  ".png"); graph_path
png(graph_path, width = Width_HalfCol*1.5, height = Width_HalfCol*0.5, units = 'in', res = 300,
    pointsize = AxisTxFontSizeSize)
grid.draw(p)
dev.off()
openfile(file = graph_path)


# 5. calculate expectation of posterior probability of acceleration ====
library(treeio)

tree_tb = tree %>% as_tibble()

get_all_parrent = function(node){
  l = parent(tree_tb, node)$label
  out = l
  # print(l)
  while(length(l) > 0){
    l = parent(tree_tb, l)$label
    # print(l)
    out = c(out, l)
  }
  return(out)
}

## test
# l = parent(tree_tb, "HLmelMel1")$label; l # "HLmelMel1-HLserCan1"
# l = parent(tree_tb, l)$label; l # "HLmelMel1-HLtaeGut4"
# "HLmelMel1-HLlepAsp1"
# "HLmelMel1-ficAlb2"
# "HLmelMel1-HLparMaj1"
# 
# x = get_all_parrent("HLmelMel1")

cal_expectation = function(node, input_df){
  # node = tree_tb$label[1]
  # input_df = run1_df
  
  print(paste0('current node: ', node))
  names(input_df) = gsub("_3", "", names(input_df))
  
  # get all parent nodes (ancestry) for the current node
  all_parents = try(get_all_parrent(node))
  
  # get posterior probablity of the current node
  tmp_node = input_df %>% select(sym(node)) %>% as.matrix()
  
  # calculate the sum of posterior probability of the ancestry of the node
  # ptm = proc.time()
  tmp = input_df %>% select(all_of(all_parents)) %>% as.matrix()
  tmp_sum = apply(tmp, 1, sum)
  # proc.time()-ptm # elaspsed: 1.845  (all)
  
  # cal expectation
  tmp_expect = tmp_node - tmp_sum #310886
}


## run1
ptm = proc.time()
exp_by_nodes_run1 = sapply(tree_tb$label, cal_expectation, input_df = run1_df)
proc.time()-ptm # elapsed: 72.604 

# row: cnees; col: nodes
exp_by_nodes_m_run1 = matrix(exp_by_nodes_run1, byrow = F, ncol = length(tree_tb$label), 
                             dimnames = list(run1_df$ID, tree_tb$label))
exp_sum_run1 = apply(exp_by_nodes_m_run1, 1, sum)
exp_sum_tb_run1 = exp_sum_run1 %>% as_tibble(rownames = NA) %>% 
  rownames_to_column(var = 'id') %>% rename('exp' = 'value')

## run2
ptm = proc.time()
exp_by_nodes_run2 = sapply(tree_tb$label, cal_expectation, input_df = run2_df)
proc.time()-ptm # elapsed: 72.604 (sapply); 74.169 (lapply)

# row: cnees; col: nodes
exp_by_nodes_m_run2 = matrix(exp_by_nodes_run2, byrow = F, ncol = length(tree_tb$label), 
                             dimnames = list(run2_df$ID, tree_tb$label))
exp_sum_run2 = apply(exp_by_nodes_m_run2, 1, sum)
exp_sum_tb_run2 = exp_sum_run2 %>% as_tibble(rownames = NA) %>% 
  rownames_to_column(var = 'id') %>% rename('exp' = 'value')

# length(exp_by_nodes_m_run2) # 30801031
# 30801031/89 # 346079
# length(intersect(transMatrix$id_run2, exp_sum_tb_run2$id)) # 339074
# length(setdiff(transMatrix$id_run2, run2_df$ID)) # 17097

## hist ====
nbins = 50

# exp run1
n = sum(transMatrix3$exp_run1>0, na.rm = T); n # 25520
p_exp1 = ggplot(transMatrix3, aes(x = exp_run1)) +
  geom_histogram(bins = nbins) +
  scale_y_continuous(expand = c(0,0,0.05,0), limits = c(NA, 40000)) +
  geom_vline(xintercept = 0, color = 'red') +
  annotate("text", x = 0, y = 20000,
           label = "x = 1",
           size = FontSize*1.5, hjust = 2, color = 'red') +
  ggtitle(paste0("run1\n", "# cnees (exp > 0): ", n)) +
  theme_m; p_exp1


# exp run2
n = sum(transMatrix3$exp_run2>0, na.rm = T); n # 12054
p_exp2 = ggplot(transMatrix3, aes(x = exp_run2)) +
  geom_histogram(bins = nbins) +
  scale_y_continuous(expand = c(0,0,0.05,0), limits = c(NA, 40000)) +
  geom_vline(xintercept = 0, color = 'red') +
  annotate("text", x = 0, y = 20000,
           label = "x = 1",
           size = FontSize*1.5, hjust = 2, color = 'red') +
  ggtitle(paste0("run2\n", "# cnees (exp > 0): ", n)) +
  theme_m; p_exp2


## compare run1 and run2
# pp
transMatrix3 = transMatrix %>%
  left_join(exp_sum_tb_run1, by = c("id_run1" = "id")) %>%
  left_join(exp_sum_tb_run2, by = c("id_run2" = "id"), suffix = c("_run1", "_run2")) %>%
  left_join(phylop) %>%
  left_join(run2_df2[, 2:15], by = c("id_run2" = "ID")) %>%
  mutate(logBF3_run2 = loglik_Full - loglik_Null)

# expectation of acceleration
corv = signif(cor(transMatrix3$exp_run1, transMatrix3$exp_run2, "complete.obs"), 2); corv # 0.81
p_exp_cor_b = ggplot(transMatrix3, aes(x = exp_run1, y = exp_run2)) +
  geom_point(size = 0.1, alpha = 0.5) + 
  theme_m; p_exp_cor_b

p_exp_cor = p_exp_cor_b + 
  geom_abline(intercept = 0, slope = 1, color = 'grey50', size = 0.2, linetype = 2) +
  geom_abline(intercept = 0, slope = corv, color = 'red') +
  annotate('text', x = -600, y = -200, 
           label = paste0("cor: ", corv),
           size = FontSize, hjust = 0, color = 'red'); p_exp_cor

p_exp_cor_postive = p_exp_cor_b + 
  geom_abline(intercept = 0, slope = 1, color = 'grey50', size = 0.2, linetype = 2) +
  geom_abline(intercept = 0, slope = corv, color = 'red') +
  annotate('text', x = 5, y = 11, 
           label = paste0("cor: ", corv),
           size = FontSize, hjust = 0, color = 'red') +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)); p_exp_cor_postive

p = arrangeGrob(grobs = list(p_exp1, p_exp2, p_exp_cor, p_exp_cor_postive), nrow = 2)

graph_path = paste0(project_path, "comparing_runs_exp_pp_",  Sys.Date(),  ".png"); graph_path
png(graph_path, width = Width_HalfCol*1, height = Width_HalfCol*0.5*2, units = 'in', res = 300,
    pointsize = AxisTxFontSizeSize)
grid.draw(p)
dev.off()
openfile(file = graph_path)

## compare with phyloP ====

# exp vs pval
p_exp2_pvalraw_cor = ggplot(transMatrix3, aes(x = -log10(pval_run2), y = exp_run2)) +
  geom_point(size = 0.2, alpha = 0.5) + 
  # geom_abline(intercept = 0, slope = 1, color = 'red') +
  theme_m; p_exp2_pvalraw_cor

# p_exp2_pvalraw_cor_zoomin = p_exp2_pvalraw_cor + coord_cartesian(xlim = c(NA, 0.005)); p_exp2_pvalraw_cor_zoomin

p_exp2_pvalraw_cor_zoominy = p_exp2_pvalraw_cor + coord_cartesian(ylim = c(0, NA)); p_exp2_pvalraw_cor_zoominy


# exp vs scale
corv = signif(cor(transMatrix3$scale_run2, transMatrix3$exp_run2, "complete.obs"), 2); corv # -0.54
p_exp2_scale_cor = ggplot(transMatrix3, aes(x = scale_run2, y = exp_run2)) +
  geom_point(size = 0.2, alpha = 0.5) + 
  # geom_abline(intercept = 0, slope = -1, color = 'grey50', size = 0.2, linetype = 2) +
  # geom_abline(intercept = 0, slope = corv, color = 'red') +
  annotate('text', x = 0.25, y = -700,
           label = paste0("cor: ", corv),
           size = FontSize, hjust = 0, color = 'red') +
  theme_m; p_exp2_scale_cor

# exp vs logBF3
corv = signif(cor(transMatrix3$logBF3_run2, transMatrix3$exp_run2, "complete.obs"), 2); corv # -0.54
p_exp2_bf3_cor = ggplot(transMatrix3, aes(x = logBF3_run2, y = exp_run2)) +
  geom_point(size = 0.2, alpha = 0.5) + 
  geom_vline(xintercept = 1, color = 'red', size = 0.2) +
  geom_hline(yintercept = 0, color = 'red', size = 0.2) +
  # geom_abline(intercept = 0, slope = 1, color = 'grey50', size = 0.2, linetype = 2) +
  # geom_abline(intercept = 0, slope = corv, color = 'red') +
  # annotate('text', x = 100, y = 50, 
  #          label = paste0("cor: ", corv),
  #          size = FontSize, hjust = 0, color = 'red') +
  theme_m ; p_exp2_bf3_cor

p_exp2_bf3_cor_zoomin = p_exp2_bf3_cor + coord_cartesian(ylim = c(0, NA)); p_exp2_bf3_cor_zoomin

p = arrangeGrob(grobs = list(p_exp2_pvalraw_cor, p_exp2_scale_cor, p_exp2_bf3_cor,
                             p_exp2_pvalraw_cor_zoominy, p_exp2_bf3_cor_zoomin), nrow = 2)

graph_path = paste0(project_path, "comparing_run2_phyloacc_phyloP_exp_pp_",  Sys.Date(),  ".png"); graph_path
png(graph_path, width = Width_HalfCol*1.5, height = Width_HalfCol*0.5*2, units = 'in', res = 300,
    pointsize = AxisTxFontSizeSize)
grid.draw(p)
dev.off()
openfile(file = graph_path)

### RUN1 ====
## hummingbird tips
# colnames(exp_by_nodes_m_run2)[grep("pha", colnames(exp_by_nodes_m_run2))]
hummingbird_tips = c("HLcalAnn5", "HLfloFus1", "HLphaSup1")

# run1
exp_by_nodes_m_run1_sel = exp_by_nodes_m_run1[, hummingbird_tips] 
# View(exp_by_nodes_m_run1_sel)
exp_sum_run1_sel = apply(exp_by_nodes_m_run1_sel, 1, sum)
exp_sum_tb_run1_sel = exp_sum_run1_sel %>% as_tibble(rownames = NA) %>% 
  rownames_to_column(var = 'id') %>% rename('exp' = 'value')

# exp run1
n = sum(exp_sum_tb_run1_sel$exp>0, na.rm = T); n # 76112
p_exp1_hum = ggplot(exp_sum_tb_run1_sel, aes(x = exp)) +
  geom_histogram(bins = nbins) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  geom_vline(xintercept = 0, color = 'red') +
  annotate("text", x = 0, y = 20000,
           label = "x = 1",
           size = FontSize*1.5, hjust = -0.5, color = 'red') +
  ggtitle(paste0("hummingbird tips (3)\n", "exp > 0: ", n)) +
  theme_m; p_exp1_hum

## honeyeater tips
colnames(exp_by_nodes_m_run2)[grep("gra", colnames(exp_by_nodes_m_run2))]
honeyeater_tips = c("HLlicPen1", "HLlicMelCas1", "HLphyNov1", "HLgraPic1")

# run1
exp_by_nodes_m_run1_sel = exp_by_nodes_m_run1[, honeyeater_tips] 
# View(exp_by_nodes_m_run1_sel)
exp_sum_run1_sel = apply(exp_by_nodes_m_run1_sel, 1, sum)
exp_sum_tb_run1_sel = exp_sum_run1_sel %>% as_tibble(rownames = NA) %>% 
  rownames_to_column(var = 'id') %>% rename('exp' = 'value')

# exp run1
n = sum(exp_sum_tb_run1_sel$exp>0, na.rm = T); n # 36118
p_exp1_honey = ggplot(exp_sum_tb_run1_sel, aes(x = exp)) +
  geom_histogram(bins = nbins) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  geom_vline(xintercept = 0, color = 'red') +
  annotate("text", x = 0, y = 20000,
           label = "x = 1",
           size = FontSize*1.5, hjust = 2, color = 'red') +
  ggtitle(paste0("honeyeater tips (4)\n", "exp > 0: ", n)) +
  theme_m; p_exp1_honey

## sunbird tips
colnames(exp_by_nodes_m_run2)[grep("lep", colnames(exp_by_nodes_m_run2))]
sunbird_tips = c("HLcinPul1", "HLlepAsp1")

# run1
exp_by_nodes_m_run1_sel = exp_by_nodes_m_run1[, sunbird_tips] 
# View(exp_by_nodes_m_run1_sel)
exp_sum_run1_sel = apply(exp_by_nodes_m_run1_sel, 1, sum)
exp_sum_tb_run1_sel = exp_sum_run1_sel %>% as_tibble(rownames = NA) %>% 
  rownames_to_column(var = 'id') %>% rename('exp' = 'value')

# exp run1
n = sum(exp_sum_tb_run1_sel$exp>0, na.rm = T); n # 19775
p_exp1_sun = ggplot(exp_sum_tb_run1_sel, aes(x = exp)) +
  geom_histogram(bins = nbins) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  geom_vline(xintercept = 0, color = 'red') +
  annotate("text", x = 0, y = 20000,
           label = "x = 1",
           size = FontSize*1.5, hjust = 2, color = 'red') +
  ggtitle(paste0("sunbird tips (2)\n", "exp > 0: ", n)) +
  theme_m; p_exp1_sun

## parrot tips
colnames(exp_by_nodes_m_run2)[grep("ama", colnames(exp_by_nodes_m_run2))]
parrot_tips = c("HLlorGal1", "HLtriMol2", "HLaraSol1", "HLamaAes1")

# run1
exp_by_nodes_m_run1_sel = exp_by_nodes_m_run1[, parrot_tips] 
# View(exp_by_nodes_m_run1_sel)
exp_sum_run1_sel = apply(exp_by_nodes_m_run1_sel, 1, sum)
exp_sum_tb_run1_sel = exp_sum_run1_sel %>% as_tibble(rownames = NA) %>% 
  rownames_to_column(var = 'id') %>% rename('exp' = 'value')

# exp run1
n = sum(exp_sum_tb_run1_sel$exp>0, na.rm = T); n # 36118
p_exp1_parrot= ggplot(exp_sum_tb_run1_sel, aes(x = exp)) +
  geom_histogram(bins = nbins) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  geom_vline(xintercept = 0, color = 'red') +
  annotate("text", x = 0, y = 20000,
           label = "x = 1",
           size = FontSize*1.5, hjust = 2, color = 'red') +
  ggtitle(paste0("parrot tips (4)\n", "exp > 0: ", n)) +
  theme_m; p_exp1_parrot

p = arrangeGrob(grobs = list(p_exp1_hum, p_exp1_honey, p_exp1_parrot, p_exp1_sun), nrow = 2)
graph_path = paste0(project_path, "run1_exp_clades_",  Sys.Date(),  ".png"); graph_path
png(graph_path, width = Width_HalfCol*1.2, height = Width_HalfCol*1, units = 'in', res = 300,
    pointsize = AxisTxFontSizeSize)
grid.draw(p)
dev.off()
openfile(file = graph_path)

### RUN2 ====
## hummingbird tips
# colnames(exp_by_nodes_m_run2)[grep("pha", colnames(exp_by_nodes_m_run2))]
hummingbird_tips = c("HLcalAnn5", "HLfloFus1", "HLphaSup1")
all_parents = unlist(sapply(hummingbird_tips, get_all_parrent))
sel = unique(c(hummingbird_tips, all_parents))

# run2
exp_by_nodes_m_run2_sel = exp_by_nodes_m_run2[, sel] 
# View(exp_by_nodes_m_run2_sel)
exp_sum_run2_sel = apply(exp_by_nodes_m_run2_sel, 1, sum)
exp_sum_tb_run2_sel = exp_sum_run2_sel %>% as_tibble(rownames = NA) %>% 
  rownames_to_column(var = 'id') %>% rename('exp' = 'value')

# exp run2
n = sum(exp_sum_tb_run2_sel$exp>0, na.rm = T); n # 76112
p_exp2_hum = ggplot(exp_sum_tb_run2_sel, aes(x = exp)) +
  geom_histogram(bins = nbins) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  geom_vline(xintercept = 0, color = 'red') +
  annotate("text", x = 0, y = 20000,
           label = "x = 1",
           size = FontSize*1.5, hjust = 2, color = 'red') +
  ggtitle(paste0("hummingbird tips + internal nodes (", length(sel), ")\n", "# cnees (exp > 0): ", n)) +
  theme_m; p_exp2_hum

## honeyeater tips
colnames(exp_by_nodes_m_run2)[grep("gra", colnames(exp_by_nodes_m_run2))]
honeyeater_tips = c("HLlicPen1", "HLlicMelCas1", "HLphyNov1", "HLgraPic1")
all_parents = unlist(sapply(honeyeater_tips, get_all_parrent))
sel = unique(c(honeyeater_tips, all_parents))

# run2
exp_by_nodes_m_run2_sel = exp_by_nodes_m_run2[, sel] 
# View(exp_by_nodes_m_run2_sel)
exp_sum_run2_sel = apply(exp_by_nodes_m_run2_sel, 1, sum)
exp_sum_tb_run2_sel = exp_sum_run2_sel %>% as_tibble(rownames = NA) %>% 
  rownames_to_column(var = 'id') %>% rename('exp' = 'value')

# exp run2
n = sum(exp_sum_tb_run2_sel$exp>0, na.rm = T); n # 76112
p_exp2_honey = ggplot(exp_sum_tb_run2_sel, aes(x = exp)) +
  geom_histogram(bins = nbins) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  geom_vline(xintercept = 0, color = 'red') +
  annotate("text", x = 0, y = 20000,
           label = "x = 1",
           size = FontSize*1.5, hjust = 2, color = 'red') +
  ggtitle(paste0("honeyeater tips + internal nodes (", length(sel), ")\n", "# cnees (exp > 0): ", n)) +
  theme_m; p_exp2_honey

## sunbird tips
colnames(exp_by_nodes_m_run2)[grep("lep", colnames(exp_by_nodes_m_run2))]
sunbird_tips = c("HLcinPul1", "HLlepAsp1")
all_parents = unlist(sapply(sunbird_tips, get_all_parrent))
sel = unique(c(sunbird_tips, all_parents))

# run2
exp_by_nodes_m_run2_sel = exp_by_nodes_m_run2[, sel] 
# View(exp_by_nodes_m_run2_sel)
exp_sum_run2_sel = apply(exp_by_nodes_m_run2_sel, 1, sum)
exp_sum_tb_run2_sel = exp_sum_run2_sel %>% as_tibble(rownames = NA) %>% 
  rownames_to_column(var = 'id') %>% rename('exp' = 'value')

# exp run2
n = sum(exp_sum_tb_run2_sel$exp>0, na.rm = T); n # 19775
p_exp2_sun = ggplot(exp_sum_tb_run2_sel, aes(x = exp)) +
  geom_histogram(bins = nbins) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  geom_vline(xintercept = 0, color = 'red') +
  annotate("text", x = 0, y = 20000,
           label = "x = 1",
           size = FontSize*1.5, hjust = 2, color = 'red') +
  ggtitle(paste0("sunbird tips + internal nodes (", length(sel), ")\n", "# cnees (exp > 0): ", n)) +
  theme_m; p_exp2_sun

## parrot tips
colnames(exp_by_nodes_m_run2)[grep("ama", colnames(exp_by_nodes_m_run2))]
parrot_tips = c("HLlorGal1", "HLtriMol2", "HLaraSol1", "HLamaAes1")
all_parents = unlist(sapply(parrot_tips, get_all_parrent))
sel = unique(c(parrot_tips, all_parents))

# run2
exp_by_nodes_m_run2_sel = exp_by_nodes_m_run2[, sel] 
# View(exp_by_nodes_m_run2_sel)
exp_sum_run2_sel = apply(exp_by_nodes_m_run2_sel, 1, sum)
exp_sum_tb_run2_sel = exp_sum_run2_sel %>% as_tibble(rownames = NA) %>% 
  rownames_to_column(var = 'id') %>% rename('exp' = 'value')

# exp run2
n = sum(exp_sum_tb_run2_sel$exp>0, na.rm = T); n # 36118
p_exp2_parrot= ggplot(exp_sum_tb_run2_sel, aes(x = exp)) +
  geom_histogram(bins = nbins) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  geom_vline(xintercept = 0, color = 'red') +
  annotate("text", x = 0, y = 20000,
           label = "x = 1",
           size = FontSize*1.5, hjust = 2, color = 'red') +
  ggtitle(paste0("parrot tips + internal nodes (", length(sel), ")\n", "# cnees (exp > 0): ", n)) +
  theme_m; p_exp2_parrot

## falcon tips
colnames(exp_by_nodes_m_run2)[grep("fal", colnames(exp_by_nodes_m_run2))]
falcon_tips = c("HLfalTin1", "falPer1", "falChe1")
all_parents = unlist(sapply(falcon_tips, get_all_parrent))
sel = unique(c(falcon_tips, all_parents))

# run2
exp_by_nodes_m_run2_sel = exp_by_nodes_m_run2[, sel] 
# View(exp_by_nodes_m_run2_sel)
exp_sum_run2_sel = apply(exp_by_nodes_m_run2_sel, 1, sum)
exp_sum_tb_run2_sel = exp_sum_run2_sel %>% as_tibble(rownames = NA) %>% 
  rownames_to_column(var = 'id') %>% rename('exp' = 'value')

# exp run2
n = sum(exp_sum_tb_run2_sel$exp>0, na.rm = T); n # 36118
p_exp2_falcon= ggplot(exp_sum_tb_run2_sel, aes(x = exp)) +
  geom_histogram(bins = nbins) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  geom_vline(xintercept = 0, color = 'red') +
  annotate("text", x = 0, y = 20000,
           label = "x = 1",
           size = FontSize*1.5, hjust = 2, color = 'red') +
  ggtitle(paste0("falcon tips + internal nodes (", length(sel), ")\n", "# cnees (exp > 0): ", n)) +
  theme_m; p_exp2_falcon

p = arrangeGrob(grobs = list(p_exp2_hum, p_exp2_honey, p_exp2_parrot, p_exp2_sun, p_exp2_falcon), nrow = 2)
graph_path = paste0(project_path, "run2_exp_clades_",  Sys.Date(),  ".png"); graph_path
png(graph_path, width = Width_HalfCol*1.5, height = Width_HalfCol*1, units = 'in', res = 300,
    pointsize = AxisTxFontSizeSize)
grid.draw(p)
dev.off()
openfile(file = graph_path)
