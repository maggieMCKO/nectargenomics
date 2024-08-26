
library(tidyverse) 
library(grid)
library(gridExtra)
library(RColorBrewer)
library(ggvenn)
library(UpSetR)
present = function(s){ifelse(!is.na(s), 1, 0)}

# local
FontSize = 2.5
AxisTxFontSizeSize = 5 # sup:6
AxisTxFontSizeSize_s = 5 # sup:6
AxisTitleFontSizeSize = 5 # sup:6
Factor_mmtoin = 0.0393701
Width_HalfCol = 85*Factor_mmtoin 
theme_m = theme(panel.background = element_blank(),
                plot.background = element_blank(), 
                plot.title = element_text(colour = "black", size = AxisTitleFontSizeSize, hjust = 0.5),
                plot.margin = margin(0.5, 2, 0.5, 0.5, "line"),
                panel.grid = element_blank(),
                panel.border = element_blank(), 
                axis.line = element_line(color = "black", size = 0.2),
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

set.seed(100)

# set up ====
setwd("/Users/maggie/ownCloud/Baldwin/Projects/Nectar/Seq_Data/globus/CNEEanalysis_2021/03_phyloacc/2024ver/nectar_nodollo/")
project_path = paste0(getwd(), "/convergence/")

out_path = paste0(project_path, "/plots/")
dir.create(out_path)

allcnee_great = read_tsv(paste0(getwd(), "/great/overlap/cnee_gene_assignment_basalPlusExtension_proteincoding_allcnees.bed"), 
                         col_names = c("chr", "start", "end", "id", "chr2", "start2", "end2", "gene")) %>% 
  mutate(gene = gsub("LOC107051846", "IRS1", gene),
         gene = gsub("LOC107055431", "NEDD4L", gene),
         gene = gsub("gga-mir-6690", "CTSB", gene),
         gene = gsub("gga-mir-12216", "B2M", gene),
         gene = gsub("gga-mir-12230", "SLC27A4", gene),
         gene = gsub("^Ii$", "CD74", gene))

# @ Loose ====
# nectar ====
## gene-lv convergence ====
path = paste0("/Users/maggie/ownCloud/Baldwin/Projects/Nectar/Seq_Data/globus/CNEEanalysis_2021/03_phyloacc/2024ver/nectar_nodollo/great/overlap/Loose/")
target_file = list.files(path, ".*_basalPlusExtension_\\w+.bed", full.names = T); target_file # great default assign rule

acccnees_100kb_withcnees0_loose = lapply(target_file, function(s){
  # s = target_file[1]
  tmp_name = basename(s)
  tmp_name = gsub("cnee_gene_assignment_basalPlusExtension_proteincoding_Loose_", "", tmp_name)
  tmp_name = gsub(".bed", "", tmp_name)
  print(paste0("current: ", tmp_name))
  x = read_tsv(s, col_names = c("chr", "start", "end", "id", "chr2", "start2", "end2", "gene")) %>% 
    dplyr::select(gene, id) %>%
    bind_cols('group' = tmp_name) 
}) %>% bind_rows()  %>% 
  # mutate(group = gsub("_2", "", group)) %>% 
  mutate(group = gsub("sunbirds_flowerpecker", "sunbirds", group)) %>% 
  mutate(gene = gsub("LOC107051846", "IRS1", gene),
         gene = gsub("LOC107055431", "NEDD4L", gene),
         gene = gsub("gga-mir-6690", "CTSB", gene),
         gene = gsub("gga-mir-12216", "B2M", gene),
         gene = gsub("gga-mir-12230", "SLC27A4", gene),
         gene = gsub("^Ii$", "CD74", gene))
head(acccnees_100kb_withcnees0_loose)

acccnees_100kb_withcnees_loose = acccnees_100kb_withcnees0_loose %>% 
  mutate(set = paste0("acc_", group), 
         val = 1,
         analysis = "acc") %>% 
  dplyr::select(set, gene, val, analysis, group) %>% 
  distinct()
dim(acccnees_100kb_withcnees_loose) # 25585    5
head(acccnees_100kb_withcnees_loose)
# set                gene     val analysis group    
# acc_2wayconvergent KLHL1      1 acc      congvergent2way
# acc_2wayconvergent PLXNA4     1 acc      congvergent2way
unique(acccnees_100kb_withcnees_loose$group)

## keep only sets of 4 cloades
acccnees_100kb_withcnees0_loose_sel = acccnees_100kb_withcnees0_loose %>%
  filter(group %in% c("honeyeaters", "hummingbirds", "parrots", "sunbirds" )) %>% 
  dplyr::rename(clade = group)

unique(acccnees_100kb_withcnees0_loose_sel$clade)

# acccnees_100kb_withcnees0_loose_selb = acccnees_100kb_withcnees0_loose_sel %>% 
#   dplyr::select(clade, gene) %>% distinct()

### num of gene per set (great) ====
nectar_clade_order = c("hummingbirds", "parrots", "honeyeaters", "sunbirds")
nectar_clade_order_union = c("union", nectar_clade_order)
n_union = n_distinct(acccnees_100kb_withcnees0_loose_sel$gene) # 2251 <- union

p1 = acccnees_100kb_withcnees0_loose_sel %>% 
  group_by(clade) %>% summarize(n = n_distinct(gene)) %>%
  bind_rows(tibble(clade = "union", n = n_union)) %>% 
  ggplot(aes(x = fct_relevel(clade, nectar_clade_order_union), y = n)) +
  geom_col(fill = 'steelblue', alpha = 0.7) + 
  geom_text(aes(label = n), vjust = 0, size = 2*2) +
  # facet_grid(.~bf) + 
  # scale_x_discrete(name = "logBFx >= n") +
  scale_y_continuous(expand = c(0,0,0.05,0), name = "gene count", limit = c(0, 8000)) + theme_m +
  theme(axis.title.x = element_blank()); p1

### num of gene per set (great) with >=2 ====
nectar_clade_order = c("hummingbirds", "parrots", "honeyeaters", "sunbirds")
nectar_clade_order_union = c("union", "congvergent_gteq2way", nectar_clade_order)
x = acccnees_100kb_withcnees0_loose %>% 
  filter(group %in% nectar_clade_order_union)

p1c = x %>% 
  group_by(group) %>% summarize(n = n_distinct(gene)) %>%
  ggplot(aes(x = fct_relevel(group, nectar_clade_order_union), y = n)) +
  geom_col(fill = 'steelblue', alpha = 0.7) + 
  geom_text(aes(label = n), vjust = 0, size = 2*2) +
  scale_y_continuous(expand = c(0,0,0.05,0), name = "gene count"
                     # , limit = c(0, 8000)
                     ) + theme_m +
  theme(axis.title.x = element_blank()); p1c

### num of convergent genes ====
genelv_conv_b = acccnees_100kb_withcnees0_loose_sel %>% 
  group_by(gene) %>% 
  mutate(n_clade = length(unique(clade))) %>% 
  dplyr::select(gene, n_clade) %>% 
  distinct()
table(genelv_conv_b$n_clade)
# 1    2    3    4 
# 1311 1404  995 1282 
# x = genelv_conv_b %>% 
#   filter(n_clade == 4)

## plot
genelv_conv_sum = genelv_conv_b %>%
  group_by(n_clade) %>%
  summarise(`Number of genes` = n()) %>%
  bind_cols(lab = c("one-way", "two-way", "three-way", "four-way")) %>%
  mutate(lab = fct_reorder(lab, n_clade, .desc = T))

col = alpha("#dc2943ff", 0.7)
p_gene = ggplot(genelv_conv_sum, aes(x = n_clade, y = `Number of genes`)) +
  geom_col(fill = col, color = 'black', linewidth = 0.1, width = 0.7) +
  geom_text(aes(label = `Number of genes`), size = FontSize*0.5, vjust = -0.5) + 
  scale_y_continuous(expand = c(0,0), limits = c(NA, 1500), breaks = c(1500))+
  scale_x_reverse() +
  theme_m + theme(axis.title.x = element_blank(), axis.text.y = element_blank()); p_gene

### upset ====
unique(acccnees_100kb_withcnees0_loose_sel$clade)
head(acccnees_100kb_withcnees0_loose_sel)

cneelv_conv_forUpset = acccnees_100kb_withcnees0_loose_sel %>% 
  dplyr::select(clade, gene) %>% distinct() %>% 
  mutate(val = 1)
dim(cneelv_conv_forUpset) #  12232     3

## to wide
cneelv_conv_wide = cneelv_conv_forUpset %>% 
  pivot_wider(id_cols = gene, names_from = clade, values_from = val, values_fn = present, values_fill = 0) %>%
  as.data.frame()
# need presence and absence for upset (both 0 and 1) and also a data.frame (tibble doesn't work)

col = alpha("#dc2943ff", 0.7) 
p = upset(cneelv_conv_wide, 
          nintersects = 11,
          sets = c("hummingbirds", "parrots", "honeyeaters", "sunbirds"), keep.order = TRUE,
          main.bar.color = col, 
          # number.colors = 'grey20', # doesnt work anymore after R update
          matrix.color = col,
          sets.bar.color = 'grey50',
          order.by = c("degree"), # "freq", 
          mb.ratio = c(0.6, 0.4),
          empty.intersections = "on", 
          mainbar.y.max = 2500,
          set_size.scale_max = 4500
); p

graph_path = paste0(out_path, "GeneLv_upset_nectar_Loose_", Sys.Date(),".pdf"); graph_path
pdf(graph_path, width = Width_HalfCol*1.325, height = Width_HalfCol*0.55, pointsize = AxisTxFontSizeSize)
p
dev.off()

## cnee-lv convergence ====
path = paste0(getwd(), "/phyloacc-out-02-07-2024.11-55-54/results/acc_cnee_sets/")
target_file = list.files(path, "Loose_\\w+.bed", full.names = T); target_file

acccnees_loose = lapply(target_file, function(s){
  # s = target_file[1]
  print(paste0("current: ", s))
  tmp_name = basename(s)
  tmp_name = gsub("Loose_", "", tmp_name)
  tmp_name = gsub(".bed", "", tmp_name)
  
  # load in accelerated CNEEs by window
  dataa <- read_tsv(s, col_names = c("chr", "start", "end", "id"), show_col_types = F) %>%
    bind_cols('group' = tmp_name) 
}) %>% bind_rows()

acccnees_loose_sel = acccnees_loose %>% 
  filter(group %in% c("honeyeaters", "hummingbirds", "parrots", "sunbirds" )) %>% 
  # mutate(group = gsub("sunbirds_flowerpecker", "sunbirds", group)) %>% 
  distinct()

### num of acc. cnee per set ====
nectar_clade_order = c("hummingbirds", "parrots", "honeyeaters", "sunbirds")
nectar_clade_order_union = c("union", nectar_clade_order)
n_union = n_distinct(acccnees_loose_sel$id) 

p1 = acccnees_loose_sel %>% 
  group_by(group) %>% tally() %>%
  bind_rows(tibble(group = "union", n = n_union)) %>% 
  ggplot(aes(x = fct_relevel(group, nectar_clade_order_union), y = n)) +
  geom_col(fill = 'steelblue') + 
  geom_text(aes(label = n), vjust = 0, size = 2*2) +
  # facet_grid(.~bf) + 
  # scale_x_discrete(name = "logBFx >= n") +
  scale_y_continuous(expand = c(0,0,0.05,0), name = "acc. cnee count", limits = c(0, 20000)) + theme_m +
  theme(axis.title.x = element_blank()); p1 

table(acccnees_loose_sel$group)
# honeyeaters hummingbirds      parrots     sunbirds 
#   4207         8186         8314         3418 

### num of convergent cnees  ====
cneelv_conv_b = acccnees_loose_sel %>% 
  dplyr::select(group, id) %>% distinct() %>% 
  group_by(id) %>% 
  mutate(n_clade = length(unique(group))) %>% 
  dplyr::select(id, n_clade) %>% 
  distinct()
table(cneelv_conv_b$n_clade)
# 1     2     3     4 
# 5991 3768 1398  450 

x = cneelv_conv_b %>% 
  filter(n_clade == 4) %>% 
  left_join(allcnee_great)
n_distinct(x$gene) # 617
cat(sort(unique(x$gene)), sep = "\n")

## plot
cneelv_conv_sum = cneelv_conv_b %>%
  group_by(n_clade) %>%
  summarise(`Number of cnees` = n()) %>%
  bind_cols(lab = c("1", "2", "3", "4")) %>%
  mutate(lab = fct_reorder(lab, n_clade, .desc = T))

col = alpha("#dc2943ff", 0.7)
p_cnee = ggplot(cneelv_conv_sum, aes(x = lab, y = `Number of cnees`)) +
  geom_col(fill = col, color = 'black', linewidth = 0.1, width = 0.7) +
  geom_text(aes(label = `Number of cnees`), size = FontSize*0.5, vjust = -0.5) + 
  scale_y_continuous(expand = c(0,0), limits = c(NA, 6500), breaks = c(6500)) +
  theme_m + theme(axis.title.x = element_blank(), axis.text.y = element_blank()); p_cnee

p = cowplot::plot_grid(p_cnee, p_gene, ncol= 2); p

graph_path = paste0(out_path, "VarLv_Convergent_phyloacc_cnees_gene_loose_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*0.5, height = Width_HalfCol*0.15, pointsize = AxisTxFontSizeSize, onefile = TRUE) # fig3 nectarspe
p
dev.off()


### upset ====
unique(acccnees_loose_sel$group)
head(acccnees_loose_sel)

cneelv_conv_forUpset = acccnees_loose_sel %>% 
  dplyr::select(group, id) %>% distinct() %>% 
  mutate(val = 1)
dim(cneelv_conv_forUpset) #  14310     3

# number of genes for each clade
tmp_clade = unique(acccnees_loose_sel$group)
sapply(tmp_clade, function(s){
  x = acccnees_loose_sel %>% filter(group == s)
  print(s)
  print(length(unique(x$id)))
})
# honeyeaters hummingbirds      parrots     sunbirds 
#     3063         4968         5441          838 

## to wide
cneelv_conv_wide = cneelv_conv_forUpset %>% 
  pivot_wider(id_cols = id, names_from = group, values_from = val, values_fn = present, values_fill = 0) %>%
  as.data.frame()
# need presence and absence for upset (both 0 and 1) and also a data.frame (tibble doesn't work)

col = alpha("#dc2943ff", 0.7) 
p = upset(cneelv_conv_wide, 
          nintersects = 11,
          sets = c("hummingbirds", "parrots", "honeyeaters", "sunbirds"), keep.order = TRUE,
          main.bar.color = col, 
          matrix.color = col,
          sets.bar.color = 'grey50',
          order.by = c("degree"), # "freq", 
          mb.ratio = c(0.6, 0.4),
          empty.intersections = "on", 
          mainbar.y.max = 2500,
          set_size.scale_max = 8000
); p

graph_path = paste0(out_path, "CneeLv_upset_nectar_Loose_", Sys.Date(),".pdf"); graph_path
pdf(graph_path, width = Width_HalfCol*1.325, height = Width_HalfCol*0.55, pointsize = AxisTxFontSizeSize)
p
dev.off()


### venn ====
tmp_clade = c("hummingbirds",  "honeyeaters", "parrots", "sunbirds")
l = lapply(tmp_clade, function(s){
  x = acccnees_loose_sel %>% filter(group == s) %>% distinct(id) %>% pull(id)
})
names(l) <- tmp_clade
str(l)

ggvenn::ggvenn(l, show_percentage = FALSE, fill_color = c('#717494ff', '#ad86aeff', '#db9982ff', '#edea93ff'))

# non-nectar ====
## gene-lv convergence ====
control_path = "/Users/maggie/ownCloud/Baldwin/Projects/Nectar/Seq_Data/globus/CNEEanalysis_2021/03_phyloacc/2024ver/control_nodollo/great/"

target_file = list.files(paste0(control_path, "/overlap/ControlLoose/"), 
                         ".*_basalPlusExtension_\\w+.bed", full.names = T); target_file # great default assign rule

acccnees_100kb_withcnees0_looseCon = lapply(target_file, function(s){
  # s = target_file[1]
  tmp_name = basename(s)
  tmp_name = gsub("cnee_gene_assignment_basalPlusExtension_proteincoding_ControlLoose_", "", tmp_name)
  tmp_name = gsub(".bed", "", tmp_name)
  print(paste0("current: ", tmp_name))
  x = read_tsv(s, col_names = c("chr", "start", "end", "id", "chr2", "start2", "end2", "gene")) %>% 
    dplyr::select(gene, id) %>%
    bind_cols('group' = tmp_name) 
}) %>% bind_rows()  %>% 
  mutate(group = gsub("sunbirds_flowerpecker", "sunbirds", group)) %>% 
  mutate(gene = gsub("LOC107051846", "IRS1", gene),
         gene = gsub("LOC107055431", "NEDD4L", gene),
         gene = gsub("gga-mir-6690", "CTSB", gene),
         gene = gsub("gga-mir-12216", "B2M", gene),
         gene = gsub("gga-mir-12230", "SLC27A4", gene),
         gene = gsub("^Ii$", "CD74", gene))
head(acccnees_100kb_withcnees0_looseCon)

acccnees_100kb_withcnees_looseCon = acccnees_100kb_withcnees0_looseCon %>% 
  mutate(set = paste0("acc_", group), 
         val = 1,
         analysis = "acc") %>% 
  dplyr::select(set, gene, val, analysis, group) %>% 
  distinct()
dim(acccnees_100kb_withcnees_looseCon) # 11322    5
head(acccnees_100kb_withcnees_looseCon)
# set                gene     val analysis group    
# acc_2wayconvergent KLHL1      1 acc      congvergent2way
# acc_2wayconvergent PLXNA4     1 acc      congvergent2way
unique(acccnees_100kb_withcnees_looseCon$group)

## keep only sets of 4 cloades
acccnees_100kb_withcnees0_looseCon_sel = acccnees_100kb_withcnees0_looseCon %>%
  filter(group %in% c("swifts", "falcons", "lyrebirds", "passerides" )) %>% 
  dplyr::rename(clade = group)

unique(acccnees_100kb_withcnees0_looseCon_sel$clade)
table(acccnees_100kb_withcnees0_looseCon_sel$clade)

##  keep only sets of 4 cloades, distinct genes (remove cnee redudency)
acccnees_100kb_withcnees0_looseCon_selb = acccnees_100kb_withcnees0_looseCon_sel %>% 
  dplyr::select(clade, gene) %>% distinct()
table(acccnees_100kb_withcnees0_looseCon_selb$clade)

### num of gene per set (great) ====
nonnectar_clade_order = c("swifts", "falcons", "lyrebirds", "passerides" )
nonnectar_clade_order_union = c("union", nonnectar_clade_order)
n_union = n_distinct(acccnees_100kb_withcnees0_looseCon_sel$gene) # 2251 <- union

p1con = acccnees_100kb_withcnees0_looseCon_sel %>% 
  group_by(clade) %>% summarize(n = n_distinct(gene)) %>%
  bind_rows(tibble(clade = "union", n = n_union)) %>% 
  ggplot(aes(x = fct_relevel(clade, nonnectar_clade_order_union), y = n)) +
  geom_col(fill = 'lightblue', alpha = 0.7) + 
  geom_text(aes(label = n), vjust = 0, size = 2*2) +
  scale_y_continuous(expand = c(0,0,0.05,0), name = "gene count", limit = c(0, 8000)) + theme_m +
  theme(axis.title.x = element_blank()); p1con

### num of gene per set (great) with >=2 ====
nonnectar_clade_order = nonnectar_clade_order = c("swifts", "falcons", "lyrebirds", "passerides" )
nonnectar_clade_order_union = c("union", "congvergent_gteq2way", nonnectar_clade_order)
x = acccnees_100kb_withcnees0_looseCon %>% 
  filter(group %in% nonnectar_clade_order_union)

p1c_con = x %>% 
  group_by(group) %>% summarize(n = n_distinct(gene)) %>%
  ggplot(aes(x = fct_relevel(group, nonnectar_clade_order_union), y = n)) +
  geom_col(fill = 'steelblue', alpha = 0.7) + 
  geom_text(aes(label = n), vjust = 0, size = 2*2) +
  scale_y_continuous(expand = c(0,0,0.05,0), name = "gene count"
                     # , limit = c(0, 8000)
  ) + theme_m +
  theme(axis.title.x = element_blank()); p1c_con

### num of convergent genes ====
genelv_conv_b = acccnees_100kb_withcnees0_looseCon_sel %>% 
  group_by(gene) %>% 
  mutate(n_clade = length(unique(clade))) %>% 
  dplyr::select(gene, n_clade) %>% 
  distinct()
table(genelv_conv_b$n_clade)
table(genelv_conv_b$n_clade)/2905
# 1    2    3    4 
# 1300  772  532  301 
# x = genelv_conv_b %>% 
#   filter(n_clade == 4)

## plot
genelv_conv_sum = genelv_conv_b %>%
  group_by(n_clade) %>%
  summarise(`Number of genes` = n()) %>%
  bind_cols(lab = c("one-way", "two-way", "three-way", "four-way")) %>%
  mutate(lab = fct_reorder(lab, n_clade, .desc = T))

col = alpha("grey", 0.7)
p_gene_con = ggplot(genelv_conv_sum, aes(x = n_clade, y = `Number of genes`)) +
  geom_col(fill = col, color = 'black', linewidth = 0.1, width = 0.7) +
  geom_text(aes(label = `Number of genes`), size = FontSize*0.5, vjust = -0.5) + # # sup: FontSize*0.8
  scale_y_continuous(expand = c(0,0), limits = c(NA, 1500), breaks = c(1500))+
  scale_x_reverse() +
  theme_m + theme(axis.title.x = element_blank(), axis.text.y = element_blank()); p_gene_con

### upset ====
unique(acccnees_100kb_withcnees0_looseCon_sel$clade)
head(acccnees_100kb_withcnees0_looseCon_sel)

cneelv_conv_forUpset = acccnees_100kb_withcnees0_looseCon_sel %>% 
  dplyr::select(clade, gene) %>% distinct() %>% 
  mutate(val = 1)
dim(cneelv_conv_forUpset) #  5644     3

## to wide
cneelv_conv_wide = cneelv_conv_forUpset %>% 
  pivot_wider(id_cols = gene, names_from = clade, values_from = val, values_fn = present, values_fill = 0) %>%
  as.data.frame()
# need presence and absence for upset (both 0 and 1) and also a data.frame (tibble doesn't work)

colcon = alpha("grey40", 1) 
p = upset(cneelv_conv_wide, 
          nintersects = 11,
          sets = c("swifts", "falcons", "lyrebirds", "passerides"), keep.order = TRUE,
          main.bar.color = colcon, 
          matrix.color = colcon,
          sets.bar.color = 'grey50',
          order.by = c("degree"), # "freq", 
          mb.ratio = c(0.6, 0.4),
          empty.intersections = "on", 
          mainbar.y.max = 2500,
          set_size.scale_max = 4500
); p

graph_path = paste0(out_path, "GeneLv_upset_nonnectar_Loose_", Sys.Date(),".pdf"); graph_path
pdf(graph_path, width = Width_HalfCol*1.325, height = Width_HalfCol*0.55, pointsize = AxisTxFontSizeSize)
p
dev.off()

## cnee-lv convergence ====
path = paste0("/Users/maggie/ownCloud/Baldwin/Projects/Nectar/Seq_Data/globus/CNEEanalysis_2021/03_phyloacc/2024ver/control_nodollo/phyloacc-out-02-20-2024.01-26-12/results/acc_cnee_sets/")
target_file = list.files(path, "ControlLoose_\\w+.bed", full.names = T); target_file

acccnees_loose_con = lapply(target_file, function(s){
  # s = target_file[8]
  print(paste0("current: ", s))
  tmp_name = basename(s)
  tmp_name = gsub("ControlLoose_", "", tmp_name)
  tmp_name = gsub(".bed", "", tmp_name)
  
  # load in accelerated CNEEs by window
  dataa <- read_tsv(s, col_names = c("chr", "start", "end", "id"), show_col_types = F) %>%
    bind_cols('group' = tmp_name) 
}) %>% bind_rows()
unique(acccnees_loose_con$group)

acccnees_loose_con_sel = acccnees_loose_con %>% 
  filter(group %in% c("falcons", "lyrebirds", "passerides", "swifts" )) %>% 
  distinct()

### num of convergent cnees [figS] ====
table(acccnees_loose_con_sel$group)
# falcons  lyrebirds passerides     swifts 
#      918        684       1119       1579 

cneelv_conv_loose_con_b = acccnees_loose_con_sel %>% 
  dplyr::select(group, id) %>% distinct() %>% 
  group_by(id) %>% 
  mutate(n_clade = length(unique(group))) %>% 
  dplyr::select(id, n_clade) %>% 
  distinct()
table(cneelv_conv_loose_con_b$n_clade)
max_nway = max(unique(cneelv_conv_loose_con_b$n_clade))

# control
# 1    2    3 
# 2996  568   56 

## plot
cneelv_conv_con_sum = cneelv_conv_loose_con_b %>%
  group_by(n_clade) %>%
  summarise(`Number of cnees` = n()) %>%
  bind_cols(lab = c(as.character(1:max_nway))) %>%
  mutate(lab = fct_reorder(lab, n_clade, .desc = T))

colcon = alpha("grey", 0.7)
p_cnee_con = ggplot(cneelv_conv_con_sum, aes(x = lab, y = `Number of cnees`)) +
  geom_col(fill = colcon, color = 'black', linewidth = 0.1, width = 0.7) +
  geom_text(aes(label = `Number of cnees`), size = FontSize*0.5, vjust = -0.5) + 
  scale_y_continuous(expand = c(0,0), limits = c(NA, 6000), breaks = c(6000)) +
  theme_m + theme(axis.title.x = element_blank(), axis.text.y = element_blank()); p_cnee_con

p = cowplot::plot_grid(p_cnee_con, p_gene_con, ncol= 2); p

graph_path = paste0(out_path, "VarLv_Convergent_phyloacc_cnees_loosecon_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*0.5, height = Width_HalfCol*0.15, pointsize = AxisTxFontSizeSize, onefile = TRUE) 
p
dev.off()


### upset ====
unique(acccnees_loose_con_sel$group)
head(acccnees_loose_con_sel)

cneelv_conv_forUpset_con = acccnees_loose_con_sel %>% 
  dplyr::select(group, id) %>% distinct() %>% 
  mutate(val = 1)
dim(cneelv_conv_forUpset_con) # 2074    3

# number of genes for each clade
tmp_clade = unique(acccnees_loose_con_sel$group)
sapply(tmp_clade, function(s){
  x = acccnees_loose_con_sel %>% filter(group == s)
  print(s)
  print(length(unique(x$id)))
})
# falcons  lyrebirds passerides     swifts 
#    1072        892       1233       1771 


cneelv_conv_con_wide = cneelv_conv_forUpset_con %>% 
  pivot_wider(id_cols = id, names_from = group, values_from = val, values_fn = present, values_fill = 0) %>%
  as.data.frame()
# need presence and absence for upset (both 0 and 1) and also a data.frame (tibble doesn't work)

colcon = alpha("grey40", 1) 
p = upset(cneelv_conv_con_wide, 
          nintersects = 11,
          sets = c("swifts", "falcons", "lyrebirds", "passerides"), keep.order = TRUE,
          main.bar.color = colcon,
          matrix.color = colcon,
          sets.bar.color = "grey50",
          # number.colors = 'grey20', # doesnt work anymore after R update
          order.by = c("degree"), # "freq", 
          mb.ratio = c(0.6, 0.4),
          empty.intersections = "on",
          mainbar.y.max = 2500,
          set_size.scale_max = 8000
); p

graph_path = paste0(out_path, "CneeLv_upset_nonnectar_ControlLoose_", Sys.Date(),".pdf"); graph_path
pdf(graph_path, width = Width_HalfCol*1.325, height = Width_HalfCol*0.55, pointsize = AxisTxFontSizeSize)
p
dev.off()



# binomial - cnee-lv ====
t = 191902 # updated 2024 nectar_nodollo
# sum(score$logbf3 >0) # 191902

# assume the probabilities of control is the background rate, what's the probability of drawing t CNEEs and observe x (nectar) cnees convergently
# 2 way
k2n = c(195,263,495,222,403,2190) # nectar
k2c =  c(177,82,39,219,113,143)    # control
# p_tmp2 = mean(k2c)/t
# 3 way
k3n =  c(228,186,346,638) # nectar
k3c = c(35,91,64,28)    # control
# p_tmp3 = mean(k3c)/t
# 4 way
k4n = 450
k4c = 39 # 0 but just say 0.5
# p_tmp4 = mean(k4c)/t

## [pointing down] layout 2 ave ====
# 2 way
p_tmp2 = mean(k2c)/t
# 3 way
p_tmp3 = mean(k3c)/t
# 4 way
p_tmp4 = mean(k4c)/t

bi_fun = function(k, ptmp){binom.test(k, t, p = ptmp, alternative = "two.sided", conf.level = 0.95)$p.value} # changed to two-sided from greater on 20231114

# 2 way
tb2 = tibble(k = 0:t, `probability` = dbinom(0:t, size = t, prob = p_tmp2), class = 'empirical background rate') %>% 
  mutate(y = probability) %>% 
  mutate(level = '2 way')

y_po = 0 # bottom
tb2n = tibble(k = k2n, y = y_po, probability = k2n/t,  class = 'nectar') %>% 
  rowwise() %>% 
  mutate(p_val = bi_fun(k, p_tmp2),
         sig = ifelse(p_val <= 0.05, 'Y', 'N')) 
tb2con = tibble(k = k2c, y = y_po, probability = k2c/t,  class = 'core non-nectar') %>% 
  mutate(sig = 'N') 
tb2c = tb2n %>% bind_rows(tb2con) %>% 
  mutate(level = '2 way')
tb2c$sig = factor(tb2c$sig, c('Y', 'N'))
tb2c$class = factor(tb2c$class, c('nectar', 'core non-nectar'))

max(tb2c$k)

pb2 = ggplot(tb2c, aes(x = k, y = y, color = class)) +
  geom_point(data = tb2, size = 0.25, color = alpha("grey50", 0.5)) +
  geom_segment(aes(xend = k, y = 0.0075, yend = y, linewidth = sig, alpha = sig), 
               arrow = arrow(ends = "last", angle = 20, length = unit(0.05 *3, "cm"))) +
  scale_x_continuous(expand = c(0,0,0,0), limits = c(0, 2300)) +
  scale_y_continuous(name = "probability", expand = c(0.02,0,0.1,0)) + # bottom, ?, top, ? , limits = c(0, y_po*1.1)
  scale_color_manual(values = c("#dc2943ff", "black")) +
  scale_linewidth_manual(values = c(0.5, 0.2), guide = 'none') +
  scale_alpha_manual(values = c(1, 0.8), guide = 'none') +
  ggtitle('2 way') + 
  theme_m +
  theme(axis.title.x = element_blank(),
        panel.border = element_rect(fill = 'transparent', linewidth = 0.2),
        legend.position = c(0.8, 0.5)); pb2

# 3 way
tb3 = tibble(k = 0:t, `probability` = dbinom(0:t, size = t, prob = p_tmp3), class = 'empirical background rate') %>% 
  mutate(y = probability)

y_po = 0 # bottom
tb3n = tibble(k = k3n, y = y_po, probability = k3n/t,  class = 'nectar') %>% 
  rowwise() %>% 
  mutate(p_val = bi_fun(k, p_tmp3),
         sig = ifelse(p_val <= 0.05, 'Y', 'N')) 
tb3con = tibble(k = k3c, y = y_po, probability = k3c/t,  class = 'core non-nectar') %>% 
  mutate(sig = 'N') 
tb3c = tb3n %>% 
  bind_rows(tb3con) %>% 
  mutate(level = '3 way')
tb3c$sig = factor(tb3c$sig, c('Y', 'N'))
tb3c$class = factor(tb3c$class, c('nectar', 'core non-nectar'))

max(tb3c$k)

pb3 = ggplot(tb3c, aes(x = k, y = y, color = class)) +
  geom_point(data = tb3, size = 0.25, color = alpha("grey50", 0.5)) +
  geom_segment(aes(xend = k, y = 0.01, yend = y, linewidth = sig, alpha = sig), 
               arrow = arrow(ends = "last", angle = 20, length = unit(0.05 *3, "cm"))) +
  scale_x_continuous(expand = c(0.05,0,0.05,0), limits = c(0, 700)) +
  scale_y_continuous(name = "probability", expand = c(0.02,0,0.1,0)) + # limits = c(0, y_max*1.1)) +
  scale_color_manual(values = c("#dc2943ff", "black")) +
  scale_linewidth_manual(values = c(0.5, 0.2)) +
  scale_alpha_manual(values = c(1, 0.8)) +
  ggtitle('3 way') +
  theme_m +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(fill = 'transparent', linewidth = 0.2),
        legend.position = 'none'); pb3

# 4 way
tb4 = tibble(k = 0:t, `probability` = dbinom(0:t, size = t, prob = p_tmp4), class = 'empirical background rate') %>% 
  mutate(y = probability)

y_po = 0 # bottom
tb4n = tibble(k = k4n, y = y_po, probability = k4n/t,  class = 'nectar') %>% 
  rowwise() %>% 
  mutate(p_val = bi_fun(k, p_tmp4),
         sig = ifelse(p_val <= 0.05, 'Y', 'N')) 
tb4con = tibble(k = k4c, y = y_po, probability = k4c/t,  class = 'core non-nectar') %>% 
  mutate(sig = 'N') 
tb4c = tb4n %>% 
  bind_rows(tb4con) %>% 
  mutate(level = '4 way')
tb4c$sig = factor(tb4c$sig, c('Y', 'N'))
tb4c$class = factor(tb4c$class, c('nectar', 'core non-nectar'))

max(tb4c$k)

pb4 = ggplot(tb4c, aes(x = k, y = y, color = class)) +
  geom_point(data = tb4, size = 0.25, color = alpha("grey50", 0.5)) +
  geom_segment(aes(xend = k, y = 0.015, yend = y, linewidth = sig, alpha = sig), 
               arrow = arrow(ends = "last", angle = 20, length = unit(0.05 *3, "cm"))) + # length is for arrow head
  scale_x_continuous(expand = c(0.05,0,0.05,0), limits = c(0, 500), n.breaks = 4) +
  scale_y_continuous(name = "probability", expand = c(0.02,0,0.1,0)) + # limits = c(0, y_max*1.1)) +
  scale_color_manual(values = c("#dc2943ff", "black", alpha("grey50", 0.25))) +
  scale_linewidth_manual(values = c(0.5, 0.2)) +
  scale_alpha_manual(values = c(1, 0.8)) +
  ggtitle('4 way') + 
  theme_m +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(fill = 'transparent', linewidth = 0.2),
        legend.position = 'none'); pb4

p = cowplot::plot_grid(pb2, pb3, pb4, nrow = 1, rel_widths = c(3.5,2,1.5)); p

graph_path = paste0(out_path, "binomial_CNEEs_pointingdw_layout2_Loose_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*1.5, height = Width_HalfCol*0.3, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p
dev.off()


# binomial - gene-lv (not done yet) ====
t = 191902 # all genes associated with cnees

# assume the probabilities of control is the background rate, what's the probability of drawing t CNEEs and observe x (nectar) cnees convergently
# 2 way
k2n = c(195,263,495,222,403,2190) # nectar
k2c =  c(177,82,39,219,113,143)    # control
# p_tmp2 = mean(k2c)/t
# 3 way
k3n =  c(228,186,346,638) # nectar
k3c = c(35,91,64,28)    # control
# p_tmp3 = mean(k3c)/t
# 4 way
k4n = 450
k4c = 39 # 0 but just say 0.5
# p_tmp4 = mean(k4c)/t

## [pointing down] layout 2 ave ====
# 2 way
p_tmp2 = mean(k2c)/t
# 3 way
p_tmp3 = mean(k3c)/t
# 4 way
p_tmp4 = mean(k4c)/t

bi_fun = function(k, ptmp){binom.test(k, t, p = ptmp, alternative = "two.sided", conf.level = 0.95)$p.value} # changed to two-sided from greater on 20231114

# 2 way
tb2 = tibble(k = 0:t, `probability` = dbinom(0:t, size = t, prob = p_tmp2), class = 'empirical background rate') %>% 
  mutate(y = probability) %>% 
  mutate(level = '2 way')

y_po = 0 # bottom
tb2n = tibble(k = k2n, y = y_po, probability = k2n/t,  class = 'nectar') %>% 
  rowwise() %>% 
  mutate(p_val = bi_fun(k, p_tmp2),
         sig = ifelse(p_val <= 0.05, 'Y', 'N')) 
tb2con = tibble(k = k2c, y = y_po, probability = k2c/t,  class = 'core non-nectar') %>% 
  mutate(sig = 'N') 
tb2c = tb2n %>% bind_rows(tb2con) %>% 
  mutate(level = '2 way')
tb2c$sig = factor(tb2c$sig, c('Y', 'N'))
tb2c$class = factor(tb2c$class, c('nectar', 'core non-nectar'))

max(tb2c$k)

pb2 = ggplot(tb2c, aes(x = k, y = y, color = class)) +
  geom_point(data = tb2, size = 0.25, color = alpha("grey50", 0.5)) +
  geom_segment(aes(xend = k, y = 0.0075, yend = y, linewidth = sig, alpha = sig), 
               arrow = arrow(ends = "last", angle = 20, length = unit(0.05 *3, "cm"))) +
  scale_x_continuous(expand = c(0,0,0,0), limits = c(0, 2300)) +
  scale_y_continuous(name = "probability", expand = c(0.02,0,0.1,0)) + # bottom, ?, top, ? , limits = c(0, y_po*1.1)
  scale_color_manual(values = c("#dc2943ff", "black")) +
  scale_linewidth_manual(values = c(0.5, 0.2), guide = 'none') +
  scale_alpha_manual(values = c(1, 0.8), guide = 'none') +
  ggtitle('2 way') + 
  theme_m +
  theme(axis.title.x = element_blank(),
        panel.border = element_rect(fill = 'transparent', linewidth = 0.2),
        legend.position = c(0.8, 0.5)); pb2

# 3 way
tb3 = tibble(k = 0:t, `probability` = dbinom(0:t, size = t, prob = p_tmp3), class = 'empirical background rate') %>% 
  mutate(y = probability)

y_po = 0 # bottom
tb3n = tibble(k = k3n, y = y_po, probability = k3n/t,  class = 'nectar') %>% 
  rowwise() %>% 
  mutate(p_val = bi_fun(k, p_tmp3),
         sig = ifelse(p_val <= 0.05, 'Y', 'N')) 
tb3con = tibble(k = k3c, y = y_po, probability = k3c/t,  class = 'core non-nectar') %>% 
  mutate(sig = 'N') 
tb3c = tb3n %>% 
  bind_rows(tb3con) %>% 
  mutate(level = '3 way')
tb3c$sig = factor(tb3c$sig, c('Y', 'N'))
tb3c$class = factor(tb3c$class, c('nectar', 'core non-nectar'))

max(tb3c$k)

pb3 = ggplot(tb3c, aes(x = k, y = y, color = class)) +
  geom_point(data = tb3, size = 0.25, color = alpha("grey50", 0.5)) +
  geom_segment(aes(xend = k, y = 0.01, yend = y, linewidth = sig, alpha = sig), 
               arrow = arrow(ends = "last", angle = 20, length = unit(0.05 *3, "cm"))) +
  scale_x_continuous(expand = c(0.05,0,0.05,0), limits = c(0, 700)) +
  scale_y_continuous(name = "probability", expand = c(0.02,0,0.1,0)) + # limits = c(0, y_max*1.1)) +
  scale_color_manual(values = c("#dc2943ff", "black")) +
  scale_linewidth_manual(values = c(0.5, 0.2)) +
  scale_alpha_manual(values = c(1, 0.8)) +
  ggtitle('3 way') +
  theme_m +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(fill = 'transparent', linewidth = 0.2),
        legend.position = 'none'); pb3

# 4 way
tb4 = tibble(k = 0:t, `probability` = dbinom(0:t, size = t, prob = p_tmp4), class = 'empirical background rate') %>% 
  mutate(y = probability)

y_po = 0 # bottom
tb4n = tibble(k = k4n, y = y_po, probability = k4n/t,  class = 'nectar') %>% 
  rowwise() %>% 
  mutate(p_val = bi_fun(k, p_tmp4),
         sig = ifelse(p_val <= 0.05, 'Y', 'N')) 
tb4con = tibble(k = k4c, y = y_po, probability = k4c/t,  class = 'core non-nectar') %>% 
  mutate(sig = 'N') 
tb4c = tb4n %>% 
  bind_rows(tb4con) %>% 
  mutate(level = '4 way')
tb4c$sig = factor(tb4c$sig, c('Y', 'N'))
tb4c$class = factor(tb4c$class, c('nectar', 'core non-nectar'))

max(tb4c$k)

pb4 = ggplot(tb4c, aes(x = k, y = y, color = class)) +
  geom_point(data = tb4, size = 0.25, color = alpha("grey50", 0.5)) +
  geom_segment(aes(xend = k, y = 0.015, yend = y, linewidth = sig, alpha = sig), 
               arrow = arrow(ends = "last", angle = 20, length = unit(0.05 *3, "cm"))) + # length is for arrow head
  scale_x_continuous(expand = c(0.05,0,0.05,0), limits = c(0, 500), n.breaks = 4) +
  scale_y_continuous(name = "probability", expand = c(0.02,0,0.1,0)) + # limits = c(0, y_max*1.1)) +
  scale_color_manual(values = c("#dc2943ff", "black", alpha("grey50", 0.25))) +
  scale_linewidth_manual(values = c(0.5, 0.2)) +
  scale_alpha_manual(values = c(1, 0.8)) +
  ggtitle('4 way') + 
  theme_m +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(fill = 'transparent', linewidth = 0.2),
        legend.position = 'none'); pb4

p = cowplot::plot_grid(pb2, pb3, pb4, nrow = 1, rel_widths = c(3.5,2,1.5)); p

graph_path = paste0(out_path, "binomial_CNEEs_pointingdw_layout2_Loose_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*1.5, height = Width_HalfCol*0.3, pointsize = AxisTxFontSizeSize, onefile = TRUE)
p
dev.off()


# Fisher - cnee-lv ====
nectar_clade_order = c("hummingbirds", "parrots", "honeyeaters", "sunbirds")
nonnectar_clade_order = c("swifts", "falcons", "lyrebirds", "passerides" )

tmp_acc_criteria = "Loose"

testFunc = function(n, i, test){
  if(n == 2){
    # 2 way
    n = 2
    nect_draw = combn(nectar_clade_order, n)
    nonnect_draw = combn(nonnectar_clade_order, n)
    # nectar
    nect1 = nect_draw[, i][1]
    nect2 = nect_draw[, i][2]
    nect1list = acccnees_loose_sel %>% 
      filter(group == nect1) %>% distinct(id)
    nect2list = acccnees_loose_sel %>% 
      filter(group == nect2) %>% distinct(id)
    l = list('clade1' = nect1list$id, 'clade2' = nect2list$id)
    v = gplots::venn(l)
    isect <- attr(v, "intersection")
    tbn = tibble(acc_criteria = tmp_acc_criteria, 
                 set = 'nectar', pair = i,
                 clade1 = nect1, clade2 = nect2, 
                 convergent = length(isect$`clade1:clade2`), 
                 notconvergent = length(isect$clade1)+length(isect$clade2))
    # nonnectar
    nonnect1 = nonnect_draw[, i][1]
    nonnect2 = nonnect_draw[, i][2]
    nonnect1list = acccnees_loose_con_sel %>% 
      filter(group == nonnect1) %>% distinct(id)
    nonnect2list = acccnees_loose_con_sel %>% 
      filter(group == nonnect2) %>% distinct(id)
    l = list('clade1' = nonnect1list$id, 'clade2' = nonnect2list$id)
    v = gplots::venn(l)
    isect <- attr(v, "intersection")
    tbnn = tibble(acc_criteria = tmp_acc_criteria, 
                  set = 'nonnectar', pair = i,
                  clade1 = nonnect1, clade2 = nonnect2, 
                  convergent = length(isect$`clade1:clade2`), 
                  notconvergent = length(isect$clade1)+length(isect$clade2))
    tbout = tbn %>% bind_rows(tbnn)
    
  }else if(n == 3){
    n = 3
    nect_draw = combn(nectar_clade_order, n)
    nonnect_draw = combn(nonnectar_clade_order, n)
    # nectar
    nect1 = nect_draw[, i][1]
    nect2 = nect_draw[, i][2]
    nect3 = nect_draw[, i][3]
    nect1list = acccnees_loose_sel %>% 
      filter(group == nect1) %>% distinct(id)
    nect2list = acccnees_loose_sel %>% 
      filter(group == nect2) %>% distinct(id)
    nect3list = acccnees_loose_sel %>% 
      filter(group == nect3) %>% distinct(id)
    l = list('clade1' = nect1list$id, 'clade2' = nect2list$id, 'clade3' = nect3list$id)
    v = gplots::venn(l)
    isect <- attr(v, "intersection")
    tbn = tibble(acc_criteria = tmp_acc_criteria, 
                 set = 'nectar', pair = i,
                 clade1 = nect1, clade2 = nect2, clade3 = nect3, 
                 convergent = length(isect$`clade1:clade2:clade3`), 
                 notconvergent = length(isect$clade1)+length(isect$clade2)+length(isect$clade3)+
                   length(isect$`clade1:clade2`)+length(isect$`clade1:clade3`)+length(isect$`clade2:clade3`))
    # nonnectar
    nonnect1 = nonnect_draw[, i][1]
    nonnect2 = nonnect_draw[, i][2]
    nonnect3 = nonnect_draw[, i][3]
    nonnect1list = acccnees_loose_con_sel %>% 
      filter(group == nonnect1) %>% distinct(id)
    nonnect2list = acccnees_loose_con_sel %>% 
      filter(group == nonnect2) %>% distinct(id)
    nonnect3list = acccnees_loose_con_sel %>% 
      filter(group == nonnect3) %>% distinct(id)
    l = list('clade1' = nonnect1list$id, 'clade2' = nonnect2list$id, 'clade3' = nonnect3list$id)
    v = gplots::venn(l)
    isect <- attr(v, "intersection")
    tbnn = tibble(acc_criteria = tmp_acc_criteria, 
                  set = 'nonnectar', pair = i,
                  clade1 = nonnect1, clade2 = nonnect2, clade3 = nonnect3, 
                  convergent = length(isect$`clade1:clade2:clade3`), 
                  notconvergent = length(isect$clade1)+length(isect$clade2)+length(isect$clade3)+
                    length(isect$`clade1:clade2`)+length(isect$`clade1:clade3`)+length(isect$`clade2:clade3`))
    tbout = tbn %>% bind_rows(tbnn)
  }else if (n == 4){
    # 4 way
    n = 4
    nect_draw = combn(nectar_clade_order, n)
    nonnect_draw = combn(nonnectar_clade_order, n)
    # nectar
    tmp = acccnees_loose_sel %>% 
      group_by(id) %>% summarise(n = n_distinct(group))
    tmp_converg = tmp %>% filter(n == 4)
    tmp_notconverg = tmp %>% filter(n < 4)
    tbn = tibble(acc_criteria = tmp_acc_criteria, 
                 set = 'nectar', pair = i,
                 convergent = n_distinct(tmp_converg$id), 
                 notconvergent = n_distinct(tmp_notconverg$id))
    # nonnectar
    tmp = acccnees_loose_con_sel %>% 
      group_by(id) %>% summarise(n = n_distinct(group))
    tmp_converg = tmp %>% filter(n == 4)
    tmp_notconverg = tmp %>% filter(n < 4)
    tbnn = tibble(acc_criteria = tmp_acc_criteria, 
                  set = 'nonnectar', pair = i,
                  # clade1 = nonnect1, clade2 = nonnect2, clade3 = nonnect3, 
                  convergent = n_distinct(tmp_converg$id), 
                  notconvergent = n_distinct(tmp_notconverg$id))
    tbout = tbn %>% bind_rows(tbnn)
  }
  
  print(i)
  print(test)
  # test = 'chi' # 'fisher' 'chiYates'
  
  
  # Create a 2-by-2 data frame
  data <- tbout[,c("convergent", "notconvergent")] %>% as.matrix()
  # colnames(data) <- c("Group1", "Group2")
  rownames(data) <- c("nectar", "nonnectar")
  
  # Run Fisher exact test
  result <- fisher.test(data, alternative = "two.sided", conf.int = TRUE)
  chi <- chisq.test(data, correct = FALSE)
  chiYates  <- chisq.test(data, correct = TRUE)
  prop_min = min(prop.table(data))
  
  if(n == 2){
    tbout2 = tbout %>% 
      mutate(`odds ratio` = result$estimate, 
             `Confidence Interval (95%)_1` = result$conf.int[1],
             `Confidence Interval (95%)_2` = result$conf.int[2],
             `P-Value` = ifelse(test == 'fisher', result$p.value, 
                                ifelse(test == 'chi', chi$p.value,
                                       ifelse(test == 'chiYates', chiYates$p.value, '?'))),
             test = test,
             nway = n, min_prop = prop_min)
  }else if(n ==3){
    tbout2 = tbout %>% 
      mutate(`odds ratio` = result$estimate, 
             `Confidence Interval (95%)_1` = result$conf.int[1],
             `Confidence Interval (95%)_2` = result$conf.int[2],
             `P-Value` = ifelse(test == 'fisher', result$p.value, 
                                ifelse(test == 'chi', chi$p.value,
                                       ifelse(test == 'chiYates', chiYates$p.value, '?'))),
             test = test,
             nway = n, min_prop = prop_min)
  }else if(n==4){
    tbout2 = tbout %>% 
      mutate(`odds ratio` = result$estimate, 
             `Confidence Interval (95%)_1` = result$conf.int[1],
             `Confidence Interval (95%)_2` = result$conf.int[2],
             `P-Value` = ifelse(test == 'fisher', result$p.value, 
                                ifelse(test == 'chi', chi$p.value,
                                       ifelse(test == 'chiYates', chiYates$p.value, '?'))),
             test = test,
             nway = n, min_prop = prop_min)
  }
}

n = rep(c(rep(2, each = 6), rep(3, each = 4), rep(4, each = 1)), 3)
i = rep(c(1:6, 1:4, 1), times = 3)
test = rep(c('fisher', 'chi', 'chiYates'), each = (6+4+1))
rbind(n,i,test)

df = mapply(testFunc,
            n = n,
            i = i,
            test = test, SIMPLIFY = FALSE)

# plot
dft = df %>% bind_rows() %>% 
  # filter(set == 'nectar') %>% 
  mutate(sig = ifelse(`P-Value`<0.05, 'Y', 'N'),
         `P-Value2` = ifelse(`P-Value`<10^-16, 10^-16, `P-Value`),
         lab = paste0("P = ", signif(`P-Value2`, digits = 1))) %>% 
  dplyr::select(test, acc_criteria, nway, set, pair, matches("clade"), convergent, notconvergent, 
                `odds ratio`, `Confidence Interval (95%)_1`, `Confidence Interval (95%)_2`,
                `P-Value`, `P-Value2`, lab, sig, min_prop)
dft$sig = factor(dft$sig, levels = c('N', 'Y'))

# reduce for plotting
dft2 = dft %>% 
  filter(set == 'nectar', test == 'chi')
  # filter(set == 'nectar', test == 'fisher')

# reduce 1
dft_reduce = dft %>% 
  filter(set == 'nectar') %>% 
  dplyr::select(test, acc_criteria, nway, pair, sig) %>% distinct() %>% 
  pivot_wider(id_cols = acc_criteria:pair, names_from = test, values_from = sig)

# reduced 2: for contigency table
dft_reduce_contingencytb = dft %>% 
  filter(test == 'fisher') %>% distinct()
path = paste0(out_path, 'CNEElv_contingency_table.csv')
# write_csv(dft_reduce_contingencytb, file = path)
# openfile(path)

# dft_chi = dft
# dft_fiser = dft

pl = lapply(2:4, function(ii){
  
  N = 5
  y_lim = ceiling(max(dft2$`Confidence Interval (95%)_2`)/N)*N
  
  tmp_dft = dft2 %>% filter(nway == ii)
  x_lim = max(tmp_dft$pair)
  
  p = ggplot(tmp_dft, aes(x = pair, color = sig)) +
    geom_hline(yintercept = 1, color = 'red', linetype = 2, linewidth = 0.25) +
    geom_point(aes( y = `odds ratio`), size = 0.75) +
    geom_errorbar(aes(ymin = `Confidence Interval (95%)_1`, 
                      ymax = `Confidence Interval (95%)_2`), width = 0) +
    geom_text(aes(label = lab, y = `Confidence Interval (95%)_2`), vjust = -0.5, size = FontSize*0.5) +
    scale_x_continuous(expand = c(0.1,0.1), breaks = 1:x_lim) +
    scale_y_continuous(expand = c(0,0,0.05,0), 
                       limits = c(0, y_lim),
                       breaks = seq(0, y_lim, by = 2)) +
    scale_color_manual(values = c('grey70', 'black'), drop = F) +
    ggtitle(paste0(ii, ' way')) + 
    # facet_grid(.~nway, space = 'free', scales = 'free') + 
    theme_m + 
    theme(axis.title.x = element_blank(),
          # axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(hjust = 1),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = 'none'
    )
  if( ii!=2 ){p = p + theme(axis.title.y = element_blank(),
                            axis.text.y = element_blank())}
  return(p)
})
p = cowplot::plot_grid(plotlist = pl, rel_widths = c(5,3,1.5), nrow = 1); p

graph_path = paste0(out_path, "Fisher_cneelv_", tmp_acc_criteria, "_chi_", Sys.Date(), ".pdf")
# pdf(graph_path, width = Width_HalfCol*1.5, height = Width_HalfCol*0.5, pointsize = AxisTxFontSizeSize, onefile = TRUE)
pdf(graph_path, width = Width_HalfCol*0.85, height = Width_HalfCol*0.25, pointsize = AxisTxFontSizeSize, onefile = TRUE) # smaller font
p
dev.off()



# Fisher - gene-lv ====
nectar_clade_order = c("hummingbirds", "parrots", "honeyeaters", "sunbirds")
nonnectar_clade_order = c("swifts", "falcons", "lyrebirds", "passerides" )

tmp_acc_criteria = "Loose"

testFunc = function(n, i, test){
  if(n == 2){
    # 2 way
    n = 2
    nect_draw = combn(nectar_clade_order, n)
    nonnect_draw = combn(nonnectar_clade_order, n)
    # nectar
    nect1 = nect_draw[, i][1]
    nect2 = nect_draw[, i][2]
    nect1list = acccnees_100kb_withcnees0_loose_sel %>% 
      filter(clade == nect1) %>% distinct(gene)
    nect2list = acccnees_100kb_withcnees0_loose_sel %>% 
      filter(clade == nect2) %>% distinct(gene)
    l = list('clade1' = nect1list$gene, 'clade2' = nect2list$gene)
    v = gplots::venn(l)
    isect <- attr(v, "intersection")
    tbn = tibble(acc_criteria = tmp_acc_criteria, 
                 set = 'nectar', pair = i,
                 clade1 = nect1, clade2 = nect2, 
                 convergent = length(isect$`clade1:clade2`), 
                 notconvergent = length(isect$clade1)+length(isect$clade2))
    # nonnectar
    nonnect1 = nonnect_draw[, i][1]
    nonnect2 = nonnect_draw[, i][2]
    nonnect1list = acccnees_100kb_withcnees0_looseCon_sel %>% 
      filter(clade == nonnect1) %>% distinct(gene)
    nonnect2list = acccnees_100kb_withcnees0_looseCon_sel %>% 
      filter(clade == nonnect2) %>% distinct(gene)
    l = list('clade1' = nonnect1list$gene, 'clade2' = nonnect2list$gene)
    v = gplots::venn(l)
    isect <- attr(v, "intersection")
    tbnn = tibble(acc_criteria = tmp_acc_criteria, 
                  set = 'nonnectar', pair = i,
                  clade1 = nonnect1, clade2 = nonnect2, 
                  convergent = length(isect$`clade1:clade2`), 
                  notconvergent = length(isect$clade1)+length(isect$clade2))
    tbout = tbn %>% bind_rows(tbnn)
    
  }else if(n == 3){
    n = 3
    nect_draw = combn(nectar_clade_order, n)
    nonnect_draw = combn(nonnectar_clade_order, n)
    # nectar
    nect1 = nect_draw[, i][1]
    nect2 = nect_draw[, i][2]
    nect3 = nect_draw[, i][3]
    nect1list = acccnees_100kb_withcnees0_loose_sel %>% 
      filter(clade == nect1) %>% distinct(gene)
    nect2list = acccnees_100kb_withcnees0_loose_sel %>% 
      filter(clade == nect2) %>% distinct(gene)
    nect3list = acccnees_100kb_withcnees0_loose_sel %>% 
      filter(clade == nect3) %>% distinct(gene)
    l = list('clade1' = nect1list$gene, 'clade2' = nect2list$gene, 'clade3' = nect3list$gene)
    v = gplots::venn(l)
    isect <- attr(v, "intersection")
    tbn = tibble(acc_criteria = tmp_acc_criteria, 
                 set = 'nectar', pair = i,
                 clade1 = nect1, clade2 = nect2, clade3 = nect3, 
                 convergent = length(isect$`clade1:clade2:clade3`), 
                 notconvergent = length(isect$clade1)+length(isect$clade2)+length(isect$clade3)+
                   length(isect$`clade1:clade2`)+length(isect$`clade1:clade3`)+length(isect$`clade2:clade3`))
    # nonnectar
    nonnect1 = nonnect_draw[, i][1]
    nonnect2 = nonnect_draw[, i][2]
    nonnect3 = nonnect_draw[, i][3]
    nonnect1list = acccnees_100kb_withcnees0_looseCon_sel %>% 
      filter(clade == nonnect1) %>% distinct(gene)
    nonnect2list = acccnees_100kb_withcnees0_looseCon_sel %>% 
      filter(clade == nonnect2) %>% distinct(gene)
    nonnect3list = acccnees_100kb_withcnees0_looseCon_sel %>% 
      filter(clade == nonnect3) %>% distinct(gene)
    l = list('clade1' = nonnect1list$gene, 'clade2' = nonnect2list$gene, 'clade3' = nonnect3list$gene)
    v = gplots::venn(l)
    isect <- attr(v, "intersection")
    tbnn = tibble(acc_criteria = tmp_acc_criteria, 
                  set = 'nonnectar', pair = i,
                  clade1 = nonnect1, clade2 = nonnect2, clade3 = nonnect3, 
                  convergent = length(isect$`clade1:clade2:clade3`), 
                  notconvergent = length(isect$clade1)+length(isect$clade2)+length(isect$clade3)+
                    length(isect$`clade1:clade2`)+length(isect$`clade1:clade3`)+length(isect$`clade2:clade3`))
    tbout = tbn %>% bind_rows(tbnn)
  }else if (n == 4){
    # 4 way
    n = 4
    nect_draw = combn(nectar_clade_order, n)
    nonnect_draw = combn(nonnectar_clade_order, n)
    # nectar
    tmp = acccnees_100kb_withcnees0_loose_sel %>% 
      group_by(gene) %>% summarise(n = n_distinct(clade))
    tmp_converg = tmp %>% filter(n == 4)
    tmp_notconverg = tmp %>% filter(n < 4)
    tbn = tibble(acc_criteria = tmp_acc_criteria, 
                 set = 'nectar', pair = i,
                 convergent = n_distinct(tmp_converg$gene), 
                 notconvergent = n_distinct(tmp_notconverg$gene))
    # nonnectar
    tmp = acccnees_100kb_withcnees0_looseCon_sel %>% 
      group_by(gene) %>% summarise(n = n_distinct(clade))
    tmp_converg = tmp %>% filter(n == 4)
    tmp_notconverg = tmp %>% filter(n < 4)
    tbnn = tibble(acc_criteria = tmp_acc_criteria, 
                  set = 'nonnectar', pair = i,
                  # clade1 = nonnect1, clade2 = nonnect2, clade3 = nonnect3, 
                  convergent = n_distinct(tmp_converg$gene), 
                  notconvergent = n_distinct(tmp_notconverg$gene))
    tbout = tbn %>% bind_rows(tbnn)
  }
  
  print(i)
  print(test)
  # test = 'chi' # 'fisher' 'chiYates'
  
  
  # Create a 2-by-2 data frame
  data <- tbout[,c("convergent", "notconvergent")] %>% as.matrix()
  # colnames(data) <- c("Group1", "Group2")
  rownames(data) <- c("nectar", "nonnectar")
  
  # Run Fisher exact test
  result <- fisher.test(data, alternative = "two.sided", conf.int = TRUE)
  chi <- chisq.test(data, correct = FALSE)
  chiYates  <- chisq.test(data, correct = TRUE)
  prop_min = min(prop.table(data))
  
  if(n == 2){
    tbout2 = tbout %>% 
      mutate(`odds ratio` = result$estimate, 
             `Confidence Interval (95%)_1` = result$conf.int[1],
             `Confidence Interval (95%)_2` = result$conf.int[2],
             `P-Value` = ifelse(test == 'fisher', result$p.value, 
                                ifelse(test == 'chi', chi$p.value,
                                       ifelse(test == 'chiYates', chiYates$p.value, '?'))),
             test = test,
             nway = n, min_prop = prop_min)
  }else if(n ==3){
    tbout2 = tbout %>% 
      mutate(`odds ratio` = result$estimate, 
             `Confidence Interval (95%)_1` = result$conf.int[1],
             `Confidence Interval (95%)_2` = result$conf.int[2],
             `P-Value` = ifelse(test == 'fisher', result$p.value, 
                                ifelse(test == 'chi', chi$p.value,
                                       ifelse(test == 'chiYates', chiYates$p.value, '?'))),
             test = test,
             nway = n, min_prop = prop_min)
  }else if(n==4){
    tbout2 = tbout %>% 
      mutate(`odds ratio` = result$estimate, 
             `Confidence Interval (95%)_1` = result$conf.int[1],
             `Confidence Interval (95%)_2` = result$conf.int[2],
             `P-Value` = ifelse(test == 'fisher', result$p.value, 
                                ifelse(test == 'chi', chi$p.value,
                                       ifelse(test == 'chiYates', chiYates$p.value, '?'))),
             test = test,
             nway = n, min_prop = prop_min)
  }
}

n = rep(c(rep(2, each = 6), rep(3, each = 4), rep(4, each = 1)), 3)
i = rep(c(1:6, 1:4, 1), times = 3)
test = rep(c('fisher', 'chi', 'chiYates'), each = (6+4+1))
rbind(n,i,test)

df = mapply(testFunc,
            n = n,
            i = i,
            test = test, SIMPLIFY = FALSE)

# plot
dft = df %>% bind_rows() %>% 
  # filter(set == 'nectar') %>% 
  mutate(sig = ifelse(`P-Value`<0.05, 'Y', 'N'),
         `P-Value2` = ifelse(`P-Value`<10^-16, 10^-16, `P-Value`),
         lab = paste0("P = ", signif(`P-Value2`, digits = 1))) %>% 
  dplyr::select(test, acc_criteria, nway, set, pair, matches("clade"), convergent, notconvergent, 
                `odds ratio`, `Confidence Interval (95%)_1`, `Confidence Interval (95%)_2`,
                `P-Value`, `P-Value2`, lab, sig, min_prop)
dft$sig = factor(dft$sig, levels = c('N', 'Y'))

# reduce for plotting
dft2 = dft %>% 
  filter(set == 'nectar', test == 'chi')
  # filter(set == 'nectar', test == 'fisher')

# reduce 1
dft_reduce = dft %>% 
  filter(set == 'nectar') %>% 
  dplyr::select(test, acc_criteria, nway, pair, sig) %>% distinct() %>% 
  pivot_wider(id_cols = acc_criteria:pair, names_from = test, values_from = sig)

# reduced 2: for contigency table
dft_reduce_contingencytb = dft %>% 
  filter(test == 'fisher') %>% distinct()
path = paste0(out_path, 'CNEElv_contingency_table_genelv_', tmp_acc_criteria, '.csv')
# write_csv(dft_reduce_contingencytb, file = path)
# openfile(path)

# dft_chi = dft
# dft_fiser = dft

pl = lapply(2:4, function(ii){
  
  y_lim = 6
  
  tmp_dft = dft2 %>% filter(nway == ii)
  x_lim = max(tmp_dft$pair)
  
  p = ggplot(tmp_dft, aes(x = pair, color = sig)) +
    geom_hline(yintercept = 1, color = 'red', linetype = 2, linewidth = 0.25) +
    geom_point(aes( y = `odds ratio`), size = 0.75) +
    geom_errorbar(aes(ymin = `Confidence Interval (95%)_1`, 
                      ymax = `Confidence Interval (95%)_2`), width = 0) +
    geom_text(aes(label = lab, y = `Confidence Interval (95%)_2`), vjust = -0.5, size = FontSize*0.5) +
    scale_x_continuous(expand = c(0.1,0.1), breaks = 1:x_lim) +
    scale_y_continuous(expand = c(0,0,0.05,0), 
                       limits = c(0, y_lim),
                       breaks = seq(0, y_lim, by = 2)) +
    scale_color_manual(values = c('grey70', 'black'), drop = F) +
    ggtitle(paste0(ii, ' way')) + 
    theme_m + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(hjust = 1),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = 'none'
    )
  if( ii!=2 ){p = p + theme(axis.title.y = element_blank(),
                            axis.text.y = element_blank())}
  return(p)
})
p = cowplot::plot_grid(plotlist = pl, rel_widths = c(5,3,1.5), nrow = 1); p

graph_path = paste0(out_path, "Fisher_genelv_", tmp_acc_criteria, "_chi_", Sys.Date(), ".pdf")
pdf(graph_path, width = Width_HalfCol*0.85, height = Width_HalfCol*0.25, pointsize = AxisTxFontSizeSize, onefile = TRUE) # smaller font
p
dev.off()

