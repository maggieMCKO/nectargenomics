library("tidyverse")
library("gridExtra")
library("grid")

FontSize = 5
AxisTxFontSizeSize = 8
AxisTxFontSizeSize_s = 6
AxisTitleFontSizeSize = 10
Factor_mmtoin = 0.0393701
Width_HalfCol = 98*Factor_mmtoin 

theme_m = 
  theme(
    plot.background = element_blank(), 
    plot.title = element_text(hjust = 0.5, size = AxisTitleFontSizeSize),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(), 
    panel.spacing.y = unit(0.2, "lines"),
    axis.line = element_line(color = "black", size = 0.1),
    axis.title.x = element_text(colour = "black", size = AxisTitleFontSizeSize),
    axis.title.y = element_text(colour = "black", size = AxisTitleFontSizeSize), # , hjust = 1, vjust = 0.5
    axis.text.x = element_text(colour = "black", size = AxisTxFontSizeSize),
    axis.text.y = element_text(colour = "black", size = AxisTxFontSizeSize), # t,r,b,l, hjust = 1
    axis.ticks = element_line(colour = "black", size = 0.02),
    strip.text = element_text(colour = "black", size = AxisTxFontSizeSize, 
                              margin = margin(2,0,2,0, "pt")),
    legend.key = element_blank(),
    legend.key.size = unit(8, 'pt'),
    legend.text = element_text(colour = "black", size = AxisTxFontSizeSize),
    legend.title = element_text(colour = "black", size = AxisTitleFontSizeSize),
    legend.background = element_blank(),
    legend.position = "none")


project_path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis_2021/01b_cnee_prep/")

# A. before merging ====
path = paste0(project_path, "cnee_len_ana/ce.lengths")
df = read_tsv(path, col_names = c("len", "set")) 
df1 = df #%>%
  # filter(len >=50) %>% 
  # filter(len < 10000)
df2 = df1 %>%
  group_by(set) %>%
  summarise(max = max(len), min = min(len), median = median(len), q = quantile(len, p = 0.9))

unique(df2$set)

graph_path = paste0(project_path, "/CNEEsets_length_dist_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*1.5, pointsize = AxisTxFontSizeSize, onefile = TRUE) 
par(mfrow = c(2,2))
sapply(unique(df2$set), function(s){
  df_sub = df1 %>% filter(set == s)
  hist(df_sub$len, breaks=50, main=paste0(s, " CNEEs"), xlab="Length")
  total = nrow(df_sub)
  num.gt.50 = sum(df_sub$len>50)
  num.gt.100 = sum(df_sub$len>100)
  legend("topright", legend=c(paste0("Total: ", total, "\n", 
                                     "Num > 50: ", num.gt.50, 
                                     " (", format(num.gt.50/total*100, digits = 2),  "%)", "\n", 
                                     "Num > 100: ", num.gt.100, 
                                     " (", format(num.gt.100/total*100, digits = 2),  "%)")), lwd=0, col="white", bty="n")
})
dev.off()


# A.2 before merging, trim large ones ====
df2_sub = df1 %>% filter(len <= 1000)

x = gtools::mixedsort(unique(df2_sub$set)); x
Ncol = ceiling(length(x)/2)

graph_path = paste0(project_path, "/CNEEsets_length_dist_1000_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*1.5, pointsize = AxisTxFontSizeSize, onefile = TRUE) 
par(mfrow = c(2,Ncol))
sapply(x, function(s){
  df_sub = df2_sub %>% filter(set == s)
  hist(df_sub$len, breaks=50, main=paste0(s, " CNEEs"), xlab="Length")
  total = nrow(df_sub)
  num.gt.50 = sum(df_sub$len>50)
  num.gt.100 = sum(df_sub$len>100)
  legend("topright", 
         legend=c(paste0("Total: ", total, "\n", 
                         "Num > 50: ", num.gt.50, 
                         " (", format(num.gt.50/total*100, digits = 2),  "%)", "\n", 
                         "Num > 100: ", num.gt.100, 
                         " (", format(num.gt.100/total*100, digits = 2),  "%)")), 
         lwd=0, col="white", bty="n")
})
dev.off()


# B.1. (no filter) after merging ====
tmp_files = list.files(paste0(project_path, "cnee_len_ana"), "galGal6_all_merged_\\d+bp.lengths", full.names = T); tmp_files
cneeLen = lapply(tmp_files, function(s){
  read_tsv(s, col_names = "len") %>% bind_cols("ori" = basename(s))
}) %>% bind_rows() %>% mutate(ori = gsub(".*_(\\d+bp).lengths", "\\1", ori))

x = gtools::mixedsort(unique(cneeLen$ori)); x
Ncol = ceiling(length(x)/3); Ncol

graph_path = paste0(project_path, "/CNEEsets_length_dist_afterMerging_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*1.5, pointsize = AxisTxFontSizeSize, onefile = TRUE) 
par(mfrow = c(3, Ncol))
sapply(x, function(s){
  df_sub = cneeLen %>% filter(ori == s)
  hist(df_sub$len, breaks=50, main=paste0(s, " CNEEs"), xlab="Length")
  total = nrow(df_sub)
  num.gt.50 = sum(df_sub$len>50)
  num.gt.100 = sum(df_sub$len>100)
  legend("topright", legend=c(paste0("Total: ", total, "\n", 
                                     "Num > 50: ", num.gt.50, 
                                     " (", format(num.gt.50/total*100, digits = 2),  "%)", "\n", 
                                     "Num > 100: ", num.gt.100, 
                                     " (", format(num.gt.100/total*100, digits = 2),  "%)")), lwd=0, col="white", bty="n")
})
dev.off()


# B.1.2 (no filter) remove large ones ====
cneeLen_sub = cneeLen %>% filter(len <= 1000)

x = gtools::mixedsort(unique(cneeLen_sub$ori)); x
Ncol = ceiling(length(x)/3)

graph_path = paste0(project_path, "/CNEEsets_length_dist_afterMerging_1000_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*1.5, pointsize = AxisTxFontSizeSize, onefile = TRUE) 
par(mfrow = c(3,Ncol))
sapply(x, function(s){
  df_sub = cneeLen_sub %>% filter(ori == s)
  hist(df_sub$len, breaks=50, main=paste0(s, " CNEEs"), xlab="Length")
  total = nrow(df_sub)
  num.gt.50 = sum(df_sub$len>50)
  num.gt.100 = sum(df_sub$len>100)
  legend("topright", 
         legend=c(paste0("Total: ", total, "\n", 
                         "Num > 50: ", num.gt.50, 
                         " (", format(num.gt.50/total*100, digits = 2),  "%)", "\n", 
                         "Num > 100: ", num.gt.100, 
                         " (", format(num.gt.100/total*100, digits = 2),  "%)")), 
         lwd=0, col="white", bty="n")
})
dev.off()



# B.1.3 (no filter) filtering: keep 50 to 1000 bp ====
cneeLen_sub = cneeLen %>% filter(len <= 1000 & len > 50)

x = gtools::mixedsort(unique(cneeLen_sub$ori)); x
Ncol = ceiling(length(x)/3)

graph_path = paste0(project_path, "/CNEEsets_length_dist_afterMerging_50_1000_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*1.5, pointsize = AxisTxFontSizeSize, onefile = TRUE) 
par(mfrow = c(3,Ncol))
sapply(x, function(s){
  df_sub = cneeLen_sub %>% filter(ori == s)
  hist(df_sub$len, breaks=50, main=paste0(s, " CNEEs"), xlab="Length")
  total = nrow(df_sub)
  num.gt.50 = sum(df_sub$len>50)
  num.gt.100 = sum(df_sub$len>100)
  legend("topright", 
         legend=c(paste0("Total: ", total, "\n", 
                         "Num > 50: ", num.gt.50, 
                         " (", format(num.gt.50/total*100, digits = 2),  "%)", "\n", 
                         "Num > 100: ", num.gt.100, 
                         " (", format(num.gt.100/total*100, digits = 2),  "%)")), 
         lwd=0, col="white", bty="n")
})
dev.off()


# B.1.4 (no filter) filtering: keep 100 to 1000 bp ====
cneeLen_sub = cneeLen %>% filter(len <= 1000 & len > 100)

x = gtools::mixedsort(unique(cneeLen_sub$ori)); x
Ncol = ceiling(length(x)/3)

graph_path = paste0(project_path, "/CNEEsets_length_dist_afterMerging_100_1000_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*1.5, pointsize = AxisTxFontSizeSize, onefile = TRUE) 
par(mfrow = c(3,Ncol))
sapply(x, function(s){
  df_sub = cneeLen_sub %>% filter(ori == s)
  hist(df_sub$len, breaks=50, main=paste0(s, " CNEEs"), xlab="Length")
  total = nrow(df_sub)
  num.gt.50 = sum(df_sub$len>50)
  num.gt.100 = sum(df_sub$len>100)
  legend("topright", 
         legend=c(paste0("Total: ", total, "\n", 
                         "Num > 50: ", num.gt.50, 
                         " (", format(num.gt.50/total*100, digits = 2),  "%)", "\n", 
                         "Num > 100: ", num.gt.100, 
                         " (", format(num.gt.100/total*100, digits = 2),  "%)")), 
         lwd=0, col="white", bty="n")
})
dev.off()


# B.2 (no filter) look at the total counts of cnees =====

p1 = cneeLen %>% 
  group_by(ori) %>%
  summarise(total = n()) %>%
  mutate(ori = fct_relevel(ori, paste0(0:10, "bp") ) ) %>%
  ggplot(aes(x = ori, y = total)) + 
  geom_col(color = 'black', fill = NA, size = 0.25) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  ggtitle("total cnee counts after merging") +
  coord_cartesian(ylim = c(1800000, NA)) + theme_m + theme(axis.title.x = element_blank()) ; p1

p2 = cneeLen %>% 
  filter(len <= 1000) %>% 
  group_by(ori) %>%
  summarise(total = n()) %>%
  mutate(ori = fct_relevel(ori, paste0(0:10, "bp") ) ) %>%
  ggplot(aes(x = ori, y = total)) + 
  geom_col(color = 'black', fill = NA, size = 0.25) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  ggtitle("total cnee counts after merging (len <= 1000)") +
  coord_cartesian(ylim = c(1800000, NA)) + theme_m + theme(axis.title.x = element_blank()) ; p2

p3 = cneeLen %>% 
  filter(len > 1000) %>% 
  group_by(ori) %>%
  summarise(total = n()) %>%
  mutate(ori = fct_relevel(ori, paste0(0:10, "bp") ) ) %>%
  ggplot(aes(x = ori, y = total)) + 
  geom_col(color = 'black', fill = NA, size = 0.25) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  ggtitle("removed cnee counts (len > 1000)") +
  coord_cartesian(ylim = c(3500, NA)) +
  theme_m + theme(axis.title.x = element_blank()) ; p3

low_cutoff = seq(50, 100, by = 10)
p_out = lapply(low_cutoff, function(s){
  df_sub = cneeLen %>% 
    filter(len <= 1000 & len > s) %>% 
    group_by(ori) %>%
    summarise(total = n()) %>%
    mutate(ori = fct_relevel(ori, paste0(0:10, "bp") ) ) 
  lmin = round(min(df_sub$total)) - 1000
  p3 = df_sub %>%
    ggplot(aes(x = ori, y = total)) + 
    geom_col(color = 'black', fill = NA, size = 0.25) +
    scale_y_continuous(expand = c(0,0,0.05,0)) +
    ggtitle(paste0("total cnee counts after merging\n(", s, " < len <= 1000)")) +
    coord_cartesian(ylim = c(lmin, NA)) + theme_m + theme(axis.title.x = element_blank()) ; p3
  
})

p = arrangeGrob(grobs = append (list(p1,p2, p3), p_out), ncol = 3)

graph_path = paste0(project_path, "/CNEEsets_length_dist_afterMerging_n_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*3, height = Width_HalfCol*1.5, pointsize = AxisTxFontSizeSize, onefile = TRUE) 
grid.draw(p)
dev.off()


# C.1. (filter <= 1000) after merging ====
tmp_files = list.files(paste0(project_path, "cnee_len_ana"), "galGal6_all_merged_eqless1000_\\d+bp.lengths", full.names = T); tmp_files
cneeLen = lapply(tmp_files, function(s){
  read_tsv(s, col_names = "len") %>% bind_cols("ori" = basename(s))
}) %>% bind_rows() %>% mutate(ori = gsub(".*_(\\d+bp).lengths", "\\1", ori))

x = gtools::mixedsort(unique(cneeLen$ori)); x
Ncol = ceiling(length(x)/3); Ncol

graph_path = paste0(project_path, "/CNEEsets_length_dist_filter1000_afterMerging_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*1.5, pointsize = AxisTxFontSizeSize, onefile = TRUE) 
par(mfrow = c(3,Ncol))
sapply(x, function(s){
  df_sub = cneeLen %>% filter(ori == s)
  hist(df_sub$len, breaks=50, main=paste0(s, " CNEEs"), xlab="Length")
  total = nrow(df_sub)
  num.gt.50 = sum(df_sub$len>50)
  num.gt.100 = sum(df_sub$len>100)
  legend("topright", legend=c(paste0("Total: ", total, "\n", 
                                     "Num > 50: ", num.gt.50, 
                                     " (", format(num.gt.50/total*100, digits = 2),  "%)", "\n", 
                                     "Num > 100: ", num.gt.100, 
                                     " (", format(num.gt.100/total*100, digits = 2),  "%)")), lwd=0, col="white", bty="n")
})
dev.off()


# C.1.2 (filter <= 1000) remove large ones ====
cneeLen_sub = cneeLen %>% filter(len <= 1000)

x = gtools::mixedsort(unique(cneeLen_sub$ori)); x
Ncol = ceiling(length(x)/3)

graph_path = paste0(project_path, "/CNEEsets_length_dist_filter1000_afterMerging_1000_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*1.5, pointsize = AxisTxFontSizeSize, onefile = TRUE) 
par(mfrow = c(3,Ncol))
sapply(x, function(s){
  df_sub = cneeLen_sub %>% filter(ori == s)
  hist(df_sub$len, breaks=50, main=paste0(s, " CNEEs"), xlab="Length")
  total = nrow(df_sub)
  num.gt.50 = sum(df_sub$len>50)
  num.gt.100 = sum(df_sub$len>100)
  legend("topright", 
         legend=c(paste0("Total: ", total, "\n", 
                         "Num > 50: ", num.gt.50, 
                         " (", format(num.gt.50/total*100, digits = 2),  "%)", "\n", 
                         "Num > 100: ", num.gt.100, 
                         " (", format(num.gt.100/total*100, digits = 2),  "%)")), 
         lwd=0, col="white", bty="n")
})
dev.off()



# C.1.3 (filter <= 1000) filtering: keep 50 to 1000 bp ====
cneeLen_sub = cneeLen %>% filter(len <= 1000 & len > 50)

x = gtools::mixedsort(unique(cneeLen_sub$ori)); x
Ncol = ceiling(length(x)/3)

graph_path = paste0(project_path, "/CNEEsets_length_dist_filter1000_afterMerging_50_1000_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*1.5, pointsize = AxisTxFontSizeSize, onefile = TRUE) 
par(mfrow = c(3,Ncol))
sapply(x, function(s){
  df_sub = cneeLen_sub %>% filter(ori == s)
  hist(df_sub$len, breaks=50, main=paste0(s, " CNEEs"), xlab="Length")
  total = nrow(df_sub)
  num.gt.50 = sum(df_sub$len>50)
  num.gt.100 = sum(df_sub$len>100)
  legend("topright", 
         legend=c(paste0("Total: ", total, "\n", 
                         "Num > 50: ", num.gt.50, 
                         " (", format(num.gt.50/total*100, digits = 2),  "%)", "\n", 
                         "Num > 100: ", num.gt.100, 
                         " (", format(num.gt.100/total*100, digits = 2),  "%)")), 
         lwd=0, col="white", bty="n")
})
dev.off()


# C.1.4 (filter <= 1000) filtering: keep 100 to 1000 bp ====
cneeLen_sub = cneeLen %>% filter(len <= 1000 & len > 100)

x = gtools::mixedsort(unique(cneeLen_sub$ori)); x
Ncol = ceiling(length(x)/3)

graph_path = paste0(project_path, "/CNEEsets_length_dist_filter1000_afterMerging_100_1000_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*1.5, pointsize = AxisTxFontSizeSize, onefile = TRUE) 
par(mfrow = c(3,Ncol))
sapply(x, function(s){
  df_sub = cneeLen_sub %>% filter(ori == s)
  hist(df_sub$len, breaks=50, main=paste0(s, " CNEEs"), xlab="Length")
  total = nrow(df_sub)
  num.gt.50 = sum(df_sub$len>50)
  num.gt.100 = sum(df_sub$len>100)
  legend("topright", 
         legend=c(paste0("Total: ", total, "\n", 
                         "Num > 50: ", num.gt.50, 
                         " (", format(num.gt.50/total*100, digits = 2),  "%)", "\n", 
                         "Num > 100: ", num.gt.100, 
                         " (", format(num.gt.100/total*100, digits = 2),  "%)")), 
         lwd=0, col="white", bty="n")
})
dev.off()


# C.1.5 (filter <= 1000) filtering: keep > 50  ====
cneeLen_sub = cneeLen %>% filter( len > 50)

x = gtools::mixedsort(unique(cneeLen_sub$ori)); x
Ncol = ceiling(length(x)/3)

graph_path = paste0(project_path, "/CNEEsets_length_dist_filter1000_afterMerging_50_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*1.5, pointsize = AxisTxFontSizeSize, onefile = TRUE) 
par(mfrow = c(3,Ncol))
sapply(x, function(s){
  df_sub = cneeLen_sub %>% filter(ori == s)
  hist(df_sub$len, breaks=50, main=paste0(s, " CNEEs"), xlab="Length")
  total = nrow(df_sub)
  num.gt.50 = sum(df_sub$len>50)
  num.gt.100 = sum(df_sub$len>100)
  legend("topright", 
         legend=c(paste0("Total: ", total, "\n", 
                         "Num > 50: ", num.gt.50, 
                         " (", format(num.gt.50/total*100, digits = 2),  "%)", "\n", 
                         "Num > 100: ", num.gt.100, 
                         " (", format(num.gt.100/total*100, digits = 2),  "%)")), 
         lwd=0, col="white", bty="n")
})
dev.off()


# C.2. (no filter) look at the total counts of cnees =====

p1 = cneeLen %>% 
  group_by(ori) %>%
  summarise(total = n()) %>%
  mutate(ori = fct_relevel(ori, paste0(0:10, "bp") ) ) %>%
  ggplot(aes(x = ori, y = total)) + 
  geom_col(color = 'black', fill = NA, size = 0.25) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  ggtitle("total cnee counts after merging") +
  coord_cartesian(ylim = c(1800000, NA)) + theme_m + theme(axis.title.x = element_blank()) ; p1

p2 = cneeLen %>% 
  filter(len <= 1000) %>% 
  group_by(ori) %>%
  summarise(total = n()) %>%
  mutate(ori = fct_relevel(ori, paste0(0:10, "bp") ) ) %>%
  ggplot(aes(x = ori, y = total)) + 
  geom_col(color = 'black', fill = NA, size = 0.25) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  ggtitle("total cnee counts after merging (len <= 1000)") +
  coord_cartesian(ylim = c(1800000, NA)) + theme_m + theme(axis.title.x = element_blank()) ; p2

p3 = cneeLen %>% 
  filter(len > 1000) %>% 
  group_by(ori) %>%
  summarise(total = n()) %>%
  mutate(ori = fct_relevel(ori, paste0(0:10, "bp") ) ) %>%
  ggplot(aes(x = ori, y = total)) + 
  geom_col(color = 'black', fill = NA, size = 0.25) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  ggtitle("removed cnee counts (len > 1000)") +
  coord_cartesian(ylim = c(500, NA)) +
  theme_m + theme(axis.title.x = element_blank()) ; p3

low_cutoff = seq(50, 100, by = 10)
p_out = lapply(low_cutoff, function(s){
  df_sub = cneeLen %>% 
    filter(len <= 1000 & len > s) %>% 
    group_by(ori) %>%
    summarise(total = n()) %>%
    mutate(ori = fct_relevel(ori, paste0(0:10, "bp") ) ) 
  lmin = round(min(df_sub$total)) - 1000
  p3 = df_sub %>%
    ggplot(aes(x = ori, y = total)) + 
    geom_col(color = 'black', fill = NA, size = 0.25) +
    scale_y_continuous(expand = c(0,0,0.05,0)) +
    ggtitle(paste0("total cnee counts after merging\n(", s, " < len <= 1000)")) +
    coord_cartesian(ylim = c(lmin, NA)) + theme_m + theme(axis.title.x = element_blank()) ; p3
  
})


p_out2 = lapply(low_cutoff, function(s){
  df_sub = cneeLen %>% 
    filter( len > s) %>% 
    group_by(ori) %>%
    summarise(total = n()) %>%
    mutate(ori = fct_relevel(ori, paste0(0:10, "bp") ) ) 
  lmin = round(min(df_sub$total)) - 1000
  p3 = df_sub %>%
    ggplot(aes(x = ori, y = total)) + 
    geom_col(color = 'black', fill = NA, size = 0.25) +
    scale_y_continuous(expand = c(0,0,0.05,0)) +
    ggtitle(paste0("total cnee counts after merging\n(", s, " < len)")) +
    coord_cartesian(ylim = c(lmin, NA)) + theme_m + theme(axis.title.x = element_blank()) ; p3
  
})

p = arrangeGrob(grobs = append(append(list(p1,p2, p3), p_out), p_out2), ncol = 3)

graph_path = paste0(project_path, "/CNEEsets_length_dist_filter1000_afterMerging_n_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*3, height = Width_HalfCol*1.5*2, pointsize = AxisTxFontSizeSize, onefile = TRUE) 
grid.draw(p)
dev.off()

