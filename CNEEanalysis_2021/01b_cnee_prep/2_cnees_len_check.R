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
unique(df2$set)

# all cnees
graph_path = paste0(project_path, "/CNEEsets_length_dist_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*1.5, pointsize = AxisTxFontSizeSize, onefile = TRUE) 
par(mfrow = c(2,2))
sapply(unique(df2$set), function(s){
  df_sub = df %>% filter(set == s)
  hist(df_sub$len, breaks=1000, main=paste0(s, " CNEEs"), xlab="Length")
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


# look at length > 1000 bp
graph_path = paste0(project_path, "/CNEEsets_length_dist_longer1000_", Sys.Date(), ".pdf")
t = unique(df2$set)
u = c(1000, 2500)
Ncol = 2
Nrow = ceiling(length(t)*length(u)/Ncol)
pdf(graph_path, width = Width_HalfCol*0.75*Ncol, height = Width_HalfCol*0.6*Nrow, pointsize = AxisTxFontSizeSize,
    onefile = TRUE) 
par(mfrow = c(Nrow,Ncol))
sapply(t, function(s){
  df_sub0 = df %>% filter(set == s)
  sapply(u, function(ut){
    df_sub = df_sub0 %>% filter(len > ut)
    hist(df_sub$len, breaks=1000, main=paste0(s, " CNEEs > ", ut," bp"), xlab="Length")
    abline(v = 2500, col = 'red')
    text(x=c(2500), y=5, labels=c("2500 bp"), pos=4, col="red")
    total = nrow(df_sub0)
    num.gt = nrow(df_sub)
    legend("topright", legend=c(paste0("Total: ", total, "\n", 
                                       "Num > ", ut, ": ", num.gt, 
                                       " (", format(num.gt/total*100, digits = 2),  "%)")), 
           lwd=0, col="white", bty="n")
  })
})
dev.off()

# UCSC
graph_path = paste0(project_path, "/CNEEsets_length_dist_UCSC_", Sys.Date(), ".pdf")
t = c(1000, 2000, 2500, 5000, 10000, 20000)
Ncol = 3
Nrow = ceiling(length(t)/Ncol)
pdf(graph_path, width = Width_HalfCol*0.75*Ncol, height = Width_HalfCol*0.6*Nrow, pointsize = AxisTxFontSizeSize,
    onefile = TRUE) 
par(mfrow = c(Nrow,Ncol))
sapply(t, function(s){
  df_sub0 = df %>% filter(set == "UCSC")
  df_sub = df_sub0 %>% filter(len > s)
  hist(df_sub$len, breaks=1000, main=paste0("UCSC CNEEs > ", s, " bp"), xlab="Length")
  abline(v = c(2500, 10000, 20000, 30000, 35000, 40000), col = 'red')
  text(x=c(42500), y=10, labels=c("2.5k, 10k, 20k, 30k, 35k, 40k bp"), pos=4, col="red")
  abline(v = 20000, col = 'red', lwd = 3)
  total = nrow(df_sub0)
  num.gt = nrow(df_sub)
  legend("topright", legend=c(paste0("Total: ", total, "\n", 
                                     "Num > ", s, ": ", num.gt, 
                                     " (", format(num.gt/total*100, digits = 2),  "%)")), 
         lwd=0, col="white", bty="n")
})
dev.off()



# num of cnees per set
df2 = df %>% filter((set %in% c("Craig", "Lowe", "Sackton") & len <= 2500) | 
                      (set %in% c("UCSC") & len <= 20000)) 
df3 = df2 %>% group_by(set) %>%
  summarise(`num of cnees` = n(), max = max(len), min = min(len))

p1 = df3 %>%
  ggplot(aes(x = set, y = `num of cnees`)) + 
  geom_col(color = 'black', fill = NA, size = 0.25) +
  geom_text (aes(label = `num of cnees`), vjust = -0.5, size = 3) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  ggtitle("cnee counts after filtering\n(Craig/Lowe/Sackton <=2500; UCSC <=20000)") +
  theme_m + theme(axis.title.x = element_blank()) ; p1

graph_path = paste0(project_path, "/num_cnees_filtered_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*1, height = Width_HalfCol*1, pointsize = AxisTxFontSizeSize,
    onefile = TRUE) 
p1
dev.off()


# B. filtered: 2500 and 20000 then merge ====
tmp_files = list.files(paste0(project_path, "cnee_len_ana"), "galGal6_all_filtered_merged", full.names = T); tmp_files
cneeLen = lapply(tmp_files, function(s){
  read_tsv(s, col_names = "len") %>% bind_cols("ori" = basename(s))
}) %>% bind_rows() %>% mutate(ori = gsub(".*_(\\d+bp).lengths", "\\1", ori))

# max len after merging
df4 = cneeLen %>% group_by(ori) %>%
  summarise(n = n(), max = max(len), min = min(len))

p1 = df4 %>% 
  mutate(ori = fct_relevel(ori, paste0(0:10, "bp") ) ) %>%
  ggplot(aes(x = ori, y = max)) + 
  geom_col(color = 'black', fill = NA, size = 0.25) +
  geom_text (aes(label = max), vjust = -0.5, size = 3) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  # ggtitle("max cnee lengths after filtering\n(Craig/Lowe/Sackton <=2500; UCSC <=20000)") +
  ggtitle("max cnee lengths after filtering\n(all <=2500)") +
  theme_m + theme(axis.title.x = element_blank()) ; p1

graph_path = paste0(project_path, "/max_cnee_len_filtered_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*1.5, height = Width_HalfCol*1, pointsize = AxisTxFontSizeSize,
    onefile = TRUE) 
p1
dev.off()



# look at total counts
p1 = cneeLen %>% 
  group_by(ori) %>%
  summarise(total = n()) %>%
  mutate(ori = fct_relevel(ori, paste0(0:10, "bp") ) ) %>%
  ggplot(aes(x = ori, y = total)) + 
  geom_col(color = 'black', fill = NA, size = 0.25) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  ggtitle("total cnee counts after merging") +
  coord_cartesian(ylim = c(1800000, NA)) + 
  theme_m + theme(axis.title.x = element_blank()) ; p1

p2 = cneeLen %>% 
  filter(len <= 20000) %>% 
  group_by(ori) %>%
  summarise(total = n()) %>%
  mutate(ori = fct_relevel(ori, paste0(0:10, "bp") ) ) %>%
  ggplot(aes(x = ori, y = total)) + 
  geom_col(color = 'black', fill = NA, size = 0.25) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  ggtitle("total cnee counts after merging\n(2nd filtering: len <= 20kbp)") +
  coord_cartesian(ylim = c(1800000, NA)) +
  theme_m + theme(axis.title.x = element_blank()) ; p2

p3 = cneeLen %>% 
  filter(len > 20000) %>% 
  group_by(ori) %>%
  summarise(total = n()) %>%
  mutate(ori = fct_relevel(ori, paste0(0:10, "bp") ) ) %>%
  ggplot(aes(x = ori, y = total)) + 
  geom_col(color = 'black', fill = NA, size = 0.25) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  ggtitle("total cnee counts after merging\n(2nd filtering: len > 20kbp)") +
  # coord_cartesian(ylim = c(500, NA)) +
  theme_m + theme(axis.title.x = element_blank()) ; p3

low_cutoff = seq(50, 70, by = 10)
# p_out = lapply(low_cutoff, function(s){
#   df_sub = cneeLen %>% 
#     filter(len <= 20000 & len > s) %>% 
#     group_by(ori) %>%
#     summarise(total = n()) %>%
#     mutate(ori = fct_relevel(ori, paste0(0:10, "bp") ) ) 
#   lmin = round(min(df_sub$total)) - 1000
#   p3 = df_sub %>%
#     ggplot(aes(x = ori, y = total)) + 
#     geom_col(color = 'black', fill = NA, size = 0.25) +
#     scale_y_continuous(expand = c(0,0,0.05,0)) +
#     ggtitle(paste0("total cnee counts after merging\n(", s, " < len <= 20k)")) +
#     coord_cartesian(ylim = c(lmin, NA)) + theme_m + theme(axis.title.x = element_blank()) ; p3
#   
# })

p_out2 = lapply(low_cutoff, function(s){
  df_sub = cneeLen %>% 
    filter(len > s) %>% 
    group_by(ori) %>%
    summarise(total = n()) %>%
    mutate(ori = fct_relevel(ori, paste0(0:10, "bp") ) ) 
  lmin = round(min(df_sub$total)) - 1000
  p3 = df_sub %>%
    ggplot(aes(x = ori, y = total)) + 
    geom_col(color = 'black', fill = NA, size = 0.25) +
    scale_y_continuous(expand = c(0,0,0.05,0)) +
    ggtitle(paste0("total cnee counts after merging\n(", s, " < len )")) +
    coord_cartesian(ylim = c(lmin, NA)) + theme_m + theme(axis.title.x = element_blank()) ; p3
  
})
p = arrangeGrob(grobs = append(list(p1,p2, p3), p_out2), ncol = 3)

graph_path = paste0(project_path, "/CNEEsets_length_dist_filtered_merged_n_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*3, height = Width_HalfCol*1, pointsize = AxisTxFontSizeSize, onefile = TRUE) 
grid.draw(p)
dev.off()


# look at short cnees
graph_path = paste0(project_path, "/CNEEsets_length_dist_filtered_merged5bp_short_", Sys.Date(), ".pdf")
t = c(500)
Ncol = 1
Nrow = ceiling(length(t)/Ncol)
pdf(graph_path, width = Width_HalfCol*0.75*Ncol, height = Width_HalfCol*0.75*Nrow, pointsize = AxisTxFontSizeSize,
    onefile = TRUE) 
par(mfrow = c(Nrow,Ncol))
sapply(t, function(s){
  df_sub0 = cneeLen %>% filter(ori == "5bp")
  df_sub = df_sub0 %>% filter(len <= s)
  hist(df_sub$len, breaks=seq(from = 0, to = s, by = 25), main=paste0("merged CNEEs (5bp) > ", s, " bp"), xlab="Length")
  # abline(v = c(2500, 20000, 30000, 35000, 40000), col = 'red')
  # text(x=c(42500), y=10, labels=c("2.5k, 20k, 30k, 35k, 40k bp"), pos=4, col="red")
  # abline(v = 20000, col = 'red', lwd = 3)
  total = nrow(df_sub0)
  num.gt = nrow(df_sub)
  legend("topright", legend=c(paste0("Total: ", total, "\n", 
                                     "Num > ", s, ": ", num.gt, 
                                     " (", format(num.gt/total*100, digits = 2),  "%)")), 
         lwd=0, col="white", bty="n")
})
dev.off()



# look at long cnees
graph_path = paste0(project_path, "/CNEEsets_length_dist_filtered_merged5bp_long_", Sys.Date(), ".pdf")
t = c(2500, 5000, 10000, 20000)
Ncol = 4
Nrow = ceiling(length(t)/Ncol)
pdf(graph_path, width = Width_HalfCol*0.75*Ncol, height = Width_HalfCol*0.6*Nrow, pointsize = AxisTxFontSizeSize,
    onefile = TRUE) 
par(mfrow = c(Nrow,Ncol))
sapply(t, function(s){
  df_sub0 = cneeLen %>% filter(ori == "5bp")
  df_sub = df_sub0 %>% filter(len > s)
  hist(df_sub$len, breaks=100, main=paste0("merged CNEEs (5bp) > ", s, " bp"), xlab="Length")
  abline(v = c(2500, 20000, 30000, 35000, 40000), col = 'red')
  text(x=c(42500), y=10, labels=c("2.5k, 20k, 30k, 35k, 40k bp"), pos=4, col="red")
  abline(v = 20000, col = 'red', lwd = 3)
  total = nrow(df_sub0)
  num.gt = nrow(df_sub)
  legend("topright", legend=c(paste0("Total: ", total, "\n", 
                                     "Num > ", s, ": ", num.gt, 
                                     " (", format(num.gt/total*100, digits = 2),  "%)")), 
         lwd=0, col="white", bty="n")
})
dev.off()




# len distribution
graph_path = paste0(project_path, "/CNEEsets_length_dist_filtered_merged_", Sys.Date(), ".pdf")
x = gtools::mixedsort(unique(cneeLen$ori)); x

Ncol = 3
Nrow = ceiling(length(x)/Ncol)

pdf(graph_path, width = Width_HalfCol*0.6*Ncol, height = Width_HalfCol*0.4*Nrow, pointsize = AxisTxFontSizeSize,
    onefile = TRUE) 
par(mfrow = c(Nrow,Ncol))
sapply(x, function(s){
  df_sub = cneeLen %>% filter(ori == s)
  hist(df_sub$len, breaks=100, main=paste0(s, " CNEEs"), xlab="Length")
  total = nrow(df_sub)
  num.gt.50 = sum(df_sub$len>2500)
  num.gt.100 = sum(df_sub$len>20000)
  legend("topright", legend=c(paste0("Total: ", total, "\n", 
                                     "Num > 2500: ", num.gt.50, 
                                     " (", format(num.gt.50/total*100, digits = 2),  "%)", "\n", 
                                     "Num > 20000: ", num.gt.100, 
                                     " (", format(num.gt.100/total*100, digits = 2),  "%)")), 
         lwd=0, col="white", bty="n")
})
dev.off()


