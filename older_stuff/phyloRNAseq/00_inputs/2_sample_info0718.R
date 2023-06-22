library(tidyverse)
library(grid)
library(gridExtra)
library(RColorBrewer)

input_path = paste0(getwd(), "/Seq_Data/globus/rnaseq/00_inputs/")


# 1. load data ====
## sample into
path = paste0(input_path, "ALL_info.tsv")
sample_info =  read_tsv(path)
dim(sample_info) # 191 4

sample_info_summary = sample_info %>% 
  group_by(bird, tissue) %>% tally() %>% 
  pivot_wider(id_cols = bird, names_from = tissue, values_from = n)
path = paste0(input_path, "/sample_summary.tsv")
# write_tsv(sample_info_summary, file = path, na = "")

unique(sample_info$bird)
# [1] "Anna's_hummingbird"      "zebra_finch"             "common_swift"            "New_Holland_honeyeater"  "rainbow_lorikeet"       
# [6] "cockatiel"               "purple-crowned_lorikeet"
x = sample_info  %>%  filter(bird == "purple-crowned_lorikeet")

unique(sample_info$tissue)
# [1] "liver"          "heart"          "pectoralis"     "palate"         "tongue"         "proventriculus" "duodenum"       "pancreas"       "splean"        
# [10] "kidney"         "fat"            "gonads"         "eyes"           "lung"           "intestine"      "crop"           "spleen"         "bottom_palate" 
# [19] "top_palate" 

## add batch info
# globus/rnaseq_prep/2022_oris/parse_RNAseq_all_0512.R
tmp_path = paste0(getwd(), "/Seq_Data/globus/rnaseq_prep/2022_oris/")

b1p = paste0(tmp_path, "batch1_lib_info.tsv")
b1 = read_tsv(b1p) %>% bind_cols("batch" = 1)

b2p = paste0(tmp_path, "batch2_lib_info.files.tsv")
b2 = read_tsv(b2p) %>% bind_cols("batch" = 2)

b3p = paste0(tmp_path, "batch3_lib_info.tsv")
b3 = read_tsv(b3p, col_names = c("LibraryID", "Name", "genome", "species"  )) %>% bind_cols("batch" = 3)

b4p = paste0(tmp_path, "batch4_lib_info.tsv")
b4 = read_tsv(b4p) %>% bind_cols("batch" = 4)

batches = b1 %>% bind_rows(b2) %>% bind_rows(b3) %>% bind_rows(b4) %>%
  mutate(species = ifelse(LibraryID %in% c("L81720", "L81721", "L81722"), "purple-crowned lorikeet", species ))
batches_sel = batches %>% dplyr::select(LibraryID, batch)

sample_info_update = sample_info %>% left_join(batches_sel, by = c("lib" = "LibraryID"))


# duodenum	heart	liver	pectoralis to do first and "palate"  "tongue"   "bottom_palate" "top_palate" 
sample_info_sel = sample_info_update %>% filter(tissue %in% c("duodenum",	"heart",	"liver",	"pectoralis",
                                                              "palate", "tongue", "bottom_palate", "top_palate"))

## expression (gene abundance (TPM))
path = paste0(input_path, "gene_abundance_all.tsv")
expr = read_tsv(path)
dim(expr) #  2411239 3

# take wanted tissue
expr_sub = expr %>% filter(sample_id %in% sample_info_sel$lib)
head(expr_sub)

# to wide
expra_w = expr_sub %>% 
  pivot_wider(id_cols = gene, names_from = sample_id, values_from = counts, values_fill = 0)

# 2. normalization ====
## B.2 EDGER ====
# https://f1000research.com/articles/5-1438/v2
library(edgeR)
# BiocManager::install("statmod")
tmp_matrix = expra_w[,-1] %>% as.matrix()
data_mean = tmp_matrix # row: gene; col: sample
data_mean.noNA = data_mean
data_mean.noNA[is.na(data_mean)] = 0
scale = calcNormFactors(data_mean.noNA, refColumn = 1, method="TMM")
scale = scale / scale[1]
data_mean.norm = sapply(seq(1, length(scale)), function(x) data_mean[,x] / scale[x])
data_mean.norm_log = log10(data_mean.norm+0.01)

colnames(data_mean.norm_log) = colnames(data_mean)
row.names(data_mean.norm_log) = expra_w$gene
# #print(head(data_mean.norm_log))

data_mean.norm_tb = data_mean.norm_log %>% as_tibble(rownames = 'gene')
data_mean.norm_tb_l = data_mean.norm_tb %>% 
  pivot_longer(cols = matches("^L"), names_to = 'sample', values_to = 'log10TMM') 

data_mean.norm_tb_l_anno = sample_info_sel %>% right_join(data_mean.norm_tb_l , by = c("lib" = "sample"))
unique(data_mean.norm_tb_l_anno$tissue)

path1 = paste0(input_path, "gene_abundance_TMM_selTissues_longAnno.tsv")
# write_tsv(data_mean.norm_tb_l_anno, file=path1)
# data_mean.norm_tb_l_anno = read_tsv(path1)

data_mean.norm_tb_l_anno_summary = data_mean.norm_tb_l_anno %>% 
  dplyr::select(bird, tissue) %>% distinct()


# 3. average by species and tissue ====

# wanted tissue types
expra_sub = data_mean.norm_tb_l_anno %>% 
  filter(tissue %in% c("duodenum",	"heart",	"liver",	"pectoralis", "palate", "top_palate")) %>% 
  mutate(tissue = gsub("top_palate", "palate", tissue)) # pull palate and top_palate
unique(expra_sub$bird)
unique(expra_sub$tissue)
x = expra_sub %>% dplyr::select(bird, tissue, lib) %>% distinct() %>% group_by(bird, tissue) %>% tally()

# average expression for species per tissue
expra_sub_mean = expra_sub %>% 
  group_by(bird, tissue, gene) %>% 
  summarise(mean_log10TMM = mean(log10TMM), n = n())
path1 = paste0(input_path, "/average_gene_expr_TMM_5Tissues.tsv") # 4 tissues
# write_tsv(expra_sub_mean, file = path1, na = "")
# expra_sub_mean = read_tsv(path1)

# save wide per tissue 
uniq_tissue = unique(expra_sub$tissue)
sapply(unique(expra_sub$tissue), function(s){
  # s = uniq_tissue[1]
  expra_w = expra_sub_mean %>% 
    filter(tissue == s) %>% 
    pivot_wider(id_cols = bird, names_from = gene, values_from = mean_log10TMM, values_fill = 0)
  path = paste0(input_path, "/average_gene_expr_TMM_", s,"_", Sys.Date(),".tsv") # 4 tissues
  # write_tsv(expra_w, file = path)
})
