### based on /Users/maggie/ownCloud/Baldwin/Projects/Nectar/Rscripts/phyloP_func2.R


library("tidyverse")
library("gridExtra")
library("grid")

fancy_scientific_y <- function(l) { 
  l <- format(l, scientific = TRUE) 
  l <- gsub("^(.*)e", "'\\1'e", l)
  l <- gsub("e", "%*%10^", l)
  parse(text=l) 
} 
fancy_scientific_y(80)

FontSize = 2.5
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
    plot.margin = unit(rep(1, 4), "cm"), # to increase
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

# A. LRT, CON [per chr] ====
path = paste0(project_path, "phylop_allChr_lrt_con.tsv")
df = read_tsv(path) %>%
  rename("chr" = "#chr") %>%
  filter(!grepl("#chr", chr)) # 378257 > 378215
nrow(df)
# names(df)
# unique(df$chr)

## 
nrow(df) # [galGal6_final_merged_CNEEs_named_fixchr_justchr.bed: 378230; 15 were lost during the process]
length(unique(df$name)) # 378215
range(df$pval) # 0 1
range(df$scale) # 0 1

df2 = df %>% mutate(padj = p.adjust(pval, "fdr"), 
                    pval_fix = ifelse(pval == 0, 10^-10, pval),
                    log10pval = -log10(pval_fix)) %>%
  mutate(fdr_bins = ifelse(padj > 0.1, '> 0.1', NA),
         fdr_bins = ifelse(padj <= 0.1 & padj > 0.05, '0.05-0.1', fdr_bins),
         fdr_bins = ifelse(padj <= 0.05 & padj > 0.01, '0.01-0.05', fdr_bins),
         fdr_bins = ifelse(padj <= 0.01 & padj > 0.001, '0.001 - 0.01', fdr_bins),
         fdr_bins = ifelse(padj <= 0.001, '<0.001', fdr_bins)) %>%
  mutate(len = end - start) %>% 
  # filter(len >= 50) %>% # no need
  mutate(scale_bins = ifelse(scale > 0.8, '> 0.8', NA),
         scale_bins = ifelse(scale <= 0.8 & scale > 0.6, '0.6-0.8', scale_bins),
         scale_bins = ifelse(scale <= 0.6 & scale > 0.4, '0.4-0.6', scale_bins),
         scale_bins = ifelse(scale <= 0.4 , '0-0.4', scale_bins)) %>% 
  mutate(fdr_bins = fct_relevel(fdr_bins, 
                                rev(c("> 0.1", "0.05-0.1", "0.01-0.05", "0.001 - 0.01", "<0.001" )))) %>% 
  mutate(scale_bins = fct_relevel(scale_bins, 
                                rev(c("> 0.8", '0.6-0.8', '0.4-0.6', '0-0.4'))))

df2$chr = factor(df2$chr,  gtools::mixedsort(unique(df2$chr)))


df_to_keep = df2 %>% filter(padj < 0.1) %>% select(name) # 363747
nrow(df2) - nrow(df_to_keep) # 14468 > 0k
path = paste0(project_path, "phyloP_conserved_cnees_to_keep.tsv")
write_tsv(df_to_keep, file = path, col_names = F)


hist(df_to_keep$padj, 10)

## plot ====
input = df2

# scale
x_tmp = input$scale
nbins = 20
# brks = seq(min(x_tmp), max(x_tmp),length.out=50)
res <- hist(x_tmp, breaks = nbins, plot = FALSE)
x_anno = max(x_tmp) * 0.6
# x_anno = quantile(x_tmp, 0.75)
y_anno = max(res$counts) * 0.8
p_scale = ggplot(input, aes(scale)) +
  geom_hline(yintercept = 0, color = 'gray30', size = 0.1) + 
  geom_histogram(bins = nbins, fill = 'blue') +
  scale_y_continuous(lab = fancy_scientific_y, expand = c(0,0,0.05,0), 
                     name = "num of cnees") +
  annotate("text", x = x_anno, y = y_anno, 
           label = paste0("median = ", format(median(x_tmp), digits = 2), 
                          "\nmax = ", max(x_tmp), "\nmin = ", min(x_tmp)), 
           size = FontSize, hjust = 0) +
  ggtitle('scale') +
  theme_m; p_scale

# cor: -log10(padj) vs scale
p_fdr_scale = ggplot(input, aes(y = scale, x = -log10(padj))) +
  geom_point(size = 0.5, color = "grey30", alpha = 0.5) +
  geom_hline(yintercept = 0.4, color = 'red') +
  annotate("text", y = 0.15, x = 4, 
           label = paste0("scale = 0.4"), 
           size = FontSize, hjust = 0, color = 'red') +
  scale_x_continuous(name = expression('-log'[10]*'padj')) +
  
  # geom_vline(xintercept = -log10(0.05), color = 'red') +
  # annotate("text", y = 0.1, x = 1.5, 
  #          label = paste0("pval = 0.05"), 
  #          size = FontSize, hjust = 0, color = 'red') +
  # geom_vline(xintercept = -log10(0.001), color = 'red') +
  # annotate("text", y = 0.06, x = 3.2, 
  #          label = paste0("pval = 0.001"), 
  #          size = FontSize, hjust = 0, color = 'red') +
  theme_m ; p_fdr_scale

# padj
x_tmp = input$padj
nbins = 20
hist(x_tmp, breaks = nbins)
x_anno = max(x_tmp) * 0.15
y_anno = nrow(input) * 0.4
total_n = nrow(input)
n_p_gt005 = sum(input$padj >= 0.05)
p_fdr = ggplot(input, aes(padj)) +
  geom_histogram(bins = nbins) +
  geom_hline(yintercept = 0, color = 'gray30', size = 0.1) + 
  scale_y_continuous(lab = fancy_scientific_y, expand = c(0,0,0.05,0), 
                     name = "num of cnees") +
  annotate("text", x = x_anno, y = y_anno,
           label = paste0("total: ", total_n, 
                          "\n\nNum of cnees (FDR >= 0.05):\n", n_p_gt005,
                          "(", format(n_p_gt005/total_n*100, digits = 2),  "%)"),
           size = FontSize, hjust = 0) +
  # scale_x_continuous(limits = c(50, NA)) +
  ggtitle('FDR') +
  theme_m; p_fdr

# class
df_m = input %>% mutate(
  conacc = ifelse(padj < 0.05, "padj\n<0.05", "padj\n>=0.05"),
  conacc = ifelse(padj < 0.01, "padj\n<0.01", conacc),
  conacc = ifelse(padj < 0.001, "padj\n<0.001", conacc)
) %>% group_by(conacc) %>% summarise(n = n())

# y_limt = max(df_m$n) * 1.05
pc = ggplot(df_m, aes(x = conacc, y = n)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), vjust = -0.5, size = FontSize) +
  scale_y_continuous(name = 'Count', 
                     # limits = c(0, y_limt), 
                     expand = c(0,0,0.05,0),
                     lab = fancy_scientific_y) +
  # scale_x_continuous(limits = c(50, NA)) +
  # ggtitle('padj class') +
  theme_m + theme(axis.title.x = element_blank()); pc

p_com1 = arrangeGrob(grobs = list(p_scale, p_fdr_scale, p_fdr, pc),  nrow = 2)

# plot 2 ====
# boxplot: x: fdr bin, y: scale
df_m = input %>% group_by(fdr_bins) %>% summarise(n = n()) 
pbox_fdr_scale = ggplot(input, aes(x = fdr_bins, y = scale)) +
  geom_boxplot(outlier.size = 0.2, width = 0.5, outlier.colour = "grey30") +
  geom_text(data = df_m, aes(label = n, y = 1.02), vjust = 0, size = FontSize) +
  # geom_point(size = 0.5, color = "grey30", alpha = 0.5, 
  # position = position_jitter(width = 0.2, height = 0)) + 
  scale_x_discrete(name = 'padj') +
  theme_m ; pbox_fdr_scale

# cor: -log10padj vs len
pscat_fdr_len = ggplot(input, aes(y = len, x = -log10(padj))) +
  geom_point(size = 0.5, color = "grey30", alpha = 0.5) +
  geom_hline(yintercept = 50, color = 'red') +
  annotate("text", y =30, x = 2,
           label = paste0("CNEE length = 50"),
           size = FontSize, hjust = 0, color = 'red') +
  # geom_vline(xintercept = -log10(0.05), color = 'red') +
  # annotate("text", y = 10^4, x = 1.5,
  #          label = paste0("pval = 0.05"),
  #          size = FontSize, hjust = 0, color = 'red') +
  # geom_vline(xintercept = -log10(0.001), color = 'red') +
  # annotate("text", y = 8000, x = 3.2, 
  #          label = paste0("pval = 0.001"), 
  #          size = FontSize, hjust = 0, color = 'red') +
  scale_y_log10(name = "CNEE length") +
  scale_x_continuous(name = expression('-log'[10]*'padj')) +
  theme_m ; pscat_fdr_len

# boxplot: x: score bin, y: len
df_m = input %>% group_by(scale_bins) %>% summarise(n = n())
pbox_scale_len = ggplot(input, aes(x = scale_bins, y = len)) +
  geom_boxplot(outlier.size = 0.2, width = 0.5, outlier.colour = "grey30") +
  geom_text(data = df_m, aes(label = n, y = 10000), vjust = 0, size = FontSize) +
  scale_y_log10(name = "CNEE length") +
  # geom_point(size = 0.5, color = "grey30", alpha = 0.5, 
  # position = position_jitter(width = 0.2, height = 0)) + 
  scale_x_discrete(name = 'scale') +
  theme_m ; pbox_scale_len

# boxplot: x: score bin, y: len [filtered]
df_m = input %>% filter(padj < 0.001) %>% group_by(scale_bins) %>% summarise(n = n())
sub = input  %>% filter(padj < 0.001)
pbox3 = ggplot(sub, aes(x = scale_bins, y = len)) +
  geom_boxplot(outlier.size = 0.2, width = 0.5, outlier.colour = "grey30") +
  geom_text(data = df_m, aes(label = n, y = 10000), vjust = 0, size = FontSize) +
  scale_y_log10(name = "CNEE length") +
  # geom_point(size = 0.5, color = "grey30", alpha = 0.5, 
  # position = position_jitter(width = 0.2, height = 0)) + 
  ggtitle('FDR < 0.0001') + 
  scale_x_discrete(name = 'scale') +
  theme_m ; pbox3

# cor: scale vs len
pscat_scale_fdr = ggplot(input, aes(y = len, x = scale)) +
  geom_point(size = 0.5, color = "grey30", alpha = 0.5) +
  geom_hline(yintercept = 50, color = 'red') +
  annotate("text", y = 2500, x = 0.1,
           label = paste0("CNEE length = 50"),
           size = FontSize, hjust = 0, color = 'red') +
  scale_y_continuous(name = "CNEE length") +
  # geom_vline(xintercept = -log10(0.05), color = 'red') +
  # annotate("text", y = 10^4, x = 1.5,
  #          label = paste0("pval = 0.05"),
  #          size = FontSize, hjust = 0, color = 'red') +
  # geom_vline(xintercept = -log10(0.001), color = 'red') +
  # annotate("text", y = 8000, x = 3.2, 
  #          label = paste0("pval = 0.001"), 
  #          size = FontSize, hjust = 0, color = 'red') +
  # scale_y_log10(name = "log10(CNEE length)") +
  theme_m ; pscat_scale_fdr

p_com2 = arrangeGrob(pbox_fdr_scale, pscat_fdr_len, pbox_scale_len, pscat_scale_fdr, ncol = 2) # pbox3
p = arrangeGrob(p_com1, p_com2, ncol = 1) 

graph_path = paste0(project_path, "phylop_LRT_Con_summary_", Sys.Date(),".pdf"); graph_path
pdf(graph_path, width = Width_HalfCol*2.5, height = Width_HalfCol*1.6*2, pointsize = AxisTxFontSizeSize,
    onefile = TRUE) 
grid.draw(p)
dev.off()


# per chr ====
length(unique(input$chr)) # 43

# scale
p_scale_per_chr = ggplot(input, aes(scale)) +
  geom_hline(yintercept = 0, color = 'gray30', size = 0.1) + 
  geom_histogram(bins = nbins, fill = 'blue') +
  facet_wrap(chr ~ . , ncol = 5, scales = "free") +
  scale_x_continuous(n.breaks = 5) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  # annotate("text", x = x_anno, y = y_anno, 
  #          label = paste0("median = ", median(x_tmp), "\nmax = ", 
  #                         max(x_tmp), "\nmin = ", min(x_tmp)), 
  #          size = FontSize*2, hjust = 0) +
  theme_m + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = AxisTxFontSizeSize_s),
        axis.text.y = element_text(size = AxisTxFontSizeSize_s),
        strip.text = element_text(colour = "black", 
                                  size = AxisTxFontSizeSize_s, 
                                  margin = margin(2,0,2,0, "pt"))); p_scale_per_chr

graph_path = paste0(project_path, "phylopLRT_Con_scale_histo_perChr_", Sys.Date(),".pdf"); graph_path
pdf(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*2, pointsize = AxisTxFontSizeSize,
    onefile = TRUE) 
p_scale_per_chr
dev.off()


x = input %>% filter(chr == "chrUn_NW_020110153v1")

# per chr filtered
df_m = df %>% mutate(
  conacc = ifelse(padj < 0.05, "padj\n<0.05", "padj\n>=0.05"),
  conacc = ifelse(padj < 0.01, "padj\n<0.01", conacc)
) %>% group_by(conacc, `#chr`) %>% summarise(n = n())

y_limt = max(df_m$n) * 1.05
pc2_chr = ggplot(df_m, aes(x = conacc, y = n)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), vjust = 0, size = FontSize) +
  facet_wrap(`#chr` ~ . , ncol = 8, scales = "free") +
  scale_y_continuous(name = 'Count', limits = c(0, y_limt)) +
  # scale_x_continuous(limits = c(50, NA)) +
  # ggtitle('scale <=0.4, padj < 0.05') +
  theme_m + theme(axis.title.x = element_blank(), 
                  axis.text.x = element_text(size = AxisTxFontSizeSize_s),
                  axis.text.y = element_text(size = AxisTxFontSizeSize_s),
                  strip.text = element_text(colour = "black", 
                                            size = AxisTxFontSizeSize_s, 
                                            margin = margin(2,0,2,0, "pt"))); pc2_chr

graph_path = paste0(path, "LRT_Con_padj_category_perChr_", Sys.Date(),".png"); graph_path
png(graph_path, width = Width_HalfCol*4, height = Width_HalfCol*2, units = 'in', res = 300, pointsize = AxisTxFontSizeSize)
pc2_chr
dev.off()
