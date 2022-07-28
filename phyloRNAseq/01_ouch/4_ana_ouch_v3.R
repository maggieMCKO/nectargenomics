library(tidyverse)
library(ouch)
library(parallel)


input_path = paste0(getwd(), "/Seq_Data/globus/rnaseq/00_inputs/")
project_path = paste0(getwd(), "/Seq_Data/globus/rnaseq/01_ouch/")

# input_path = paste0("/home/mpg08/mko/Nectar/analysis/rnaseq/00_inputs/")
# project_path = paste0("/home/mpg08/mko/Nectar/analysis/rnaseq/01_ouch/") 

output_path = paste0(project_path, "/output/")

path = list.files(output_path, paste0("_ouch_TMM_2022-07-19|28"), full.names = T, recursive = T); path

out_ouch_tb = lapply(path, function(s){read_tsv(s)}) %>% bind_rows()
table(out_ouch_tb$tissue)
# TMM nor
# duodenum      heart      liver     palate pectoralis 
# 14190      14147      14098      14227      14033 

# aicc ====
out_ouch_aicc = out_ouch_tb %>% 
  dplyr::select(tissue, gene, matches("aicc"), matches("pvalue"), matches("qvalue")) %>% 
  pivot_longer(cols = matches('aicc'), 
               names_pattern = "aicc_(.*)", 
               names_to = c("model"), values_to = "aicc")
out_ouch_aicc_min = out_ouch_aicc %>% 
  group_by(tissue, gene) %>% 
  slice_min(order_by = aicc, n = 1) 

# by aicc: 
table(out_ouch_aicc_min$tissue, out_ouch_aicc_min$model)
# TMM nor
#               Br    Hn
# duodenum   14179    11
# heart      14140     7
# liver      14086    12
# palate     14216    11
# pectoralis 14019    14

# by qvalue_Hg: 0
x = out_ouch_aicc_min %>% filter(qvalue_Hg<0.05); nrow(x)
# by qvalue_Hn: 322
y = out_ouch_aicc_min %>% filter(qvalue_Hn<0.05); nrow(y)
table(y$tissue)

# TMM nor
# duodenum      heart      liver     palate pectoralis 
#     76         38         76         72         60 

y_tisW = y %>% 
  pivot_wider(id_cols = gene, names_from = tissue, values_from = qvalue_Hn)


# by qvalue_HnHg: 174
z = out_ouch_aicc_min %>% filter(qvalue_HnHg<0.05); nrow(z)
table(z$tissue)
# duodenum      heart      liver     palate pectoralis 
#     36         28         40         37         33 

length(intersect(y$gene, z$gene)) # 81


### Hn vs BM ====
## BM
br_cf_Hn__br = out_ouch_tb %>% filter(qvalue_Hn >= 0.05) %>% 
  dplyr::select(tissue, gene, matches("Br")&matches("sigma|theta|qval"), matches("_Hn")&matches("sigma|theta|qval|var|alpha")) %>% 
  dplyr::rename("sigma.squared_Br" = "sigma_Br",
         "sigma.squared_Hn" = "sigma_Hn")
names(br_cf_Hn__br)
dim(br_cf_Hn__br) # 70419
table(br_cf_Hn__br$tissue)

p1s = ggplot(br_cf_Hn__br, aes(x = sigma.squared_Br, y = sigma.squared_Hn)) +
  geom_point( alpha = 0.1) + theme_classic(); p1s

p1a = ggplot(br_cf_Hn__br, aes(x = alpha_Hn)) +
  geom_histogram(bins = 10) + 
  scale_x_continuous(breaks = c(50000, 150000, 250000)) +
  theme_classic(); p1a

p1v = ggplot(br_cf_Hn__br, aes(x = var_Hn)) +
  geom_histogram(bins = 10) + theme_classic(); p1v

cor1tn = cor(br_cf_Hn__br$theta_Br, br_cf_Hn__br$theta_Hn_nectar); cor1tn #  0.784395
p1tn = ggplot(br_cf_Hn__br, aes(x = theta_Br, y = theta_Hn_nectar)) +
  geom_point( alpha = 0.1) + theme_classic(); p1tn

cor1tnn = cor(br_cf_Hn__br$theta_Br, br_cf_Hn__br$theta_Hn_nonnectar); cor1tnn #   0.9459804
p1tnn = ggplot(br_cf_Hn__br, aes(x = theta_Br, y = theta_Hn_nonnectar)) +
  geom_point( alpha = 0.1) + theme_classic(); p1tnn
# p1tnnn = ggplot(br_cf_Hn__br, aes(x = theta_Hn_nectar, y = theta_Hn_nonnectar)) +
#   geom_point( alpha = 0.1) + theme_classic(); p1tnnn


# Hn
br_cf_Hn__Hn = out_ouch_tb %>% filter(qvalue_Hn < 0.05) %>% 
  dplyr::select(tissue, gene, matches("Br")&matches("sigma|theta|qval"), matches("_Hn")&matches("sigma|theta|qval|var|alpha")) %>% 
  dplyr::rename("sigma.squared_Br" = "sigma_Br",
                "sigma.squared_Hn" = "sigma_Hn")
names(br_cf_Hn__Hn)
dim(br_cf_Hn__Hn) # 322
table(br_cf_Hn__Hn$tissue)
# duodenum      heart      liver     palate pectoralis 
#     76         38         76         72         60 


p2s = ggplot(br_cf_Hn__Hn, aes(x = sigma.squared_Br, sigma.squared_Hn)) +
  geom_point( alpha = 0.1) + theme_classic(); p2s

p2a = ggplot(br_cf_Hn__Hn, aes(x = alpha_Hn)) +
  geom_histogram(bins = 10) + theme_classic(); p2a

p2v = ggplot(br_cf_Hn__Hn, aes(x = var_Hn)) +
  geom_histogram(bins = 10) + theme_classic(); p2v

cor2tn = cor(br_cf_Hn__Hn$theta_Br, br_cf_Hn__Hn$theta_Hn_nectar); cor2tn #  0.4099429
p2tn = ggplot(br_cf_Hn__Hn, aes(x = theta_Br, y = theta_Hn_nectar)) +
  geom_point( alpha = 0.1) + theme_classic(); p2tn

cor2tn = cor(br_cf_Hn__Hn$theta_Br, br_cf_Hn__Hn$theta_Hn_nonnectar); cor2tn #  0.7429724
p2tnn = ggplot(br_cf_Hn__Hn, aes(x = theta_Br, y = theta_Hn_nonnectar)) +
  geom_point( alpha = 0.1) + theme_classic(); p2tnn

p = arrangeGrob(grobs = list(p1s, p1a, p1v, p1tn, p1tnn,
                             p2s, p2a, p2v, p2tn, p2tnn), nrow = 2)

graph_path = paste0(output_path, "/BMvsHN_diagnostics_1_", Sys.Date(), ".pdf"); 
pdf(graph_path, width = 12, height = 5)
grid.draw(p)
dev.off()
openfile(file = graph_path)


### Hn vs Hg ====
## Hg
Hg_cf_Hn__Hg = out_ouch_tb %>% filter(qvalue_HnHg >= 0.05) %>% 
  dplyr::select(tissue, gene, matches("Hg")&matches("sigma|theta|qval|var|alpha"), 
                matches("_Hn")&matches("sigma|theta|qval|var|alpha")) %>% 
  dplyr::rename("sigma.squared_Hg" = "sigma_Hg",
                "sigma.squared_Hn" = "sigma_Hn")
names(Hg_cf_Hn__Hg)
dim(Hg_cf_Hn__Hg) # 70521
table(Hg_cf_Hn__Hg$tissue)
# duodenum      heart      liver     palate pectoralis 
#    14154      14119      14058      14190      14000 

p1s = ggplot(Hg_cf_Hn__Hg, aes(x = sigma.squared_Hg, y = sigma.squared_Hn)) +
  geom_point( alpha = 0.1) + theme_classic(); p1s

p1a = ggplot(Hg_cf_Hn__Hg, aes(x = alpha_Hg, y = alpha_Hn)) +
  geom_point( alpha = 0.1) + theme_classic(); p1a

p1v = ggplot(Hg_cf_Hn__Hg, aes(x = var_Hg, y = var_Hn)) +
  geom_point( alpha = 0.1) + theme_classic(); p1v

cor1tn = cor(Hg_cf_Hn__Hg$theta_Hg, Hg_cf_Hn__Hg$theta_Hn_nectar); cor1tn # 0.7913042
p1tn = ggplot(Hg_cf_Hn__Hg, aes(x = theta_Hg, y = theta_Hn_nectar)) +
  geom_point( alpha = 0.1) + theme_classic(); p1tn
cor1tnn = cor(Hg_cf_Hn__Hg$theta_Hg, Hg_cf_Hn__Hg$theta_Hn_nonnectar); cor1tnn # 0.9532949
p1tnn = ggplot(Hg_cf_Hn__Hg, aes(x = theta_Hg, y = theta_Hn_nonnectar)) +
  geom_point( alpha = 0.1) + theme_classic(); p1tnn


# Hn
Hg_cf_Hn__Hn = out_ouch_tb %>% filter(qvalue_HnHg < 0.05) %>% 
  dplyr::select(tissue, gene, matches("Hg")&matches("sigma|theta|qval|var|alpha"), 
                matches("_Hn")&matches("sigma|theta|qval|var|alpha")) %>% 
  dplyr::rename("sigma.squared_Hg" = "sigma_Hg",
                "sigma.squared_Hn" = "sigma_Hn")
names(Hg_cf_Hn__Hn)
dim(Hg_cf_Hn__Hn) # 174
table(Hg_cf_Hn__Hn$tissue)
# duodenum      heart      liver     palate pectoralis 
#     36         28         40         37         33 

cor2s = cor(Hg_cf_Hn__Hn$sigma.squared_Hg, Hg_cf_Hn__Hn$sigma.squared_Hn); cor2s # 0.1478258
p2s = ggplot(Hg_cf_Hn__Hn, aes(x = sigma.squared_Hg, sigma.squared_Hn)) +
  geom_point( alpha = 0.1) + theme_classic(); p2s

p2a = ggplot(Hg_cf_Hn__Hn, aes(x = alpha_Hg, y = alpha_Hn)) +
  geom_point( alpha = 0.1) + theme_classic(); p2a

p2v = ggplot(Hg_cf_Hn__Hn, aes(x = var_Hg, y = var_Hn)) +
  geom_point( alpha = 0.1) + theme_classic(); p2v

cor2tn = cor(Hg_cf_Hn__Hn$theta_Hg, Hg_cf_Hn__Hn$theta_Hn_nectar); cor2tn # 0.2851806
p2tn = ggplot(Hg_cf_Hn__Hn, aes(x = theta_Hg, y = theta_Hn_nectar)) +
  geom_point( alpha = 0.1) + theme_classic(); p2tn

cor2tnn = cor(Hg_cf_Hn__Hn$theta_Hg, Hg_cf_Hn__Hn$theta_Hn_nonnectar); cor2tnn # 0.5113936
p2tnn = ggplot(Hg_cf_Hn__Hn, aes(x = theta_Hg, y = theta_Hn_nonnectar)) +
  geom_point( alpha = 0.1) + theme_classic(); p2tnn

p = arrangeGrob(grobs = list(p1s, p1a, p1v, p1tn, p1tnn,
                             p2s, p2a, p2v, p2tn, p2tnn), nrow = 2)

# dev.off()
# grid.draw(p)

graph_path = paste0(output_path, "/Hn_vs_Hg_diagnostics_1_", Sys.Date(), ".pdf"); 
pdf(graph_path, width = 12, height = 5)
grid.draw(p)
dev.off()
openfile(file = graph_path)

### anno ====
# Hg_cf_Hn__Hn
dim(Hg_cf_Hn__Hn) # 174

# transcript to gene
gene_anno = read_tsv(paste0(project_path, "../../pcoc/from_KE/isoforms.one2ones.checked.tsv"), col_names = c("transcript", "gene"))
length(unique(gene_anno$transcript)) # 31489
length(unique(gene_anno$gene)) # 9280

path2 = paste0(output_path, "Hn_vs_Hg_TMM__Hn_list_", Sys.Date(), ".tsv")
# write_tsv(Hg_cf_Hn__Hn, path2, na= "")


Hg_cf_Hn__Hn_anno = Hg_cf_Hn__Hn %>% left_join(gene_anno) %>% 
  dplyr::select(tissue, gene, transcript, matches("Hn"))
dim(Hg_cf_Hn__Hn_anno) # 310


library(biomaRt) 
ensembl = useMart("ensembl")
# datasets = listDatasets(ensembl)
# head(datasets) 
# dim(datasets) # 203 3
# datasets$dataset[grep("gal", datasets$dataset, ignore.case = T)]

ensembl_gg = useDataset("ggallus_gene_ensembl", mart=ensembl)
filters = listFilters(ensembl_gg) # query # from
filters$name[grep('symbo', filters$name)] # hgnc_symbol

attributes = listAttributes(ensembl_gg) # to
attributes$name[grep('name', attributes$name)] # wikigene_name external_gene_name
attributes$name[grep('ggallus', attributes$name)] # 

geneDescr =  getBM(attributes=c('hgnc_symbol', 'description'), 
                   filters = 'hgnc_symbol', 
                   values = Hg_cf_Hn__Hn_anno$gene, 
                   mart = ensembl_gg)
geneDescr2 = geneDescr %>%
  separate(description, into = c("description", "source"), sep = " \\[Source:", remove = F) %>%
  separate(source, into = c("source", "acc"), sep = ";Acc:", remove = T) %>%
  mutate(acc = gsub("\\]", "", acc),
         acc = gsub("HGNC:", "", acc))

go_geneDescr = Hg_cf_Hn__Hn_anno %>% 
  dplyr::left_join(geneDescr2, by = c("gene" = "hgnc_symbol")) 
nrow(go_geneDescr) # 310
good = go_geneDescr %>% filter(! is.na(description))

# look up in human
tolookagain = go_geneDescr %>% filter(is.na(description))
datasets = listDatasets(ensembl)
head(datasets)
dim(datasets) # 215 3
datasets$dataset[grep("hsa", datasets$dataset, ignore.case = T)]

ensembl_h = useDataset("hsapiens_gene_ensembl", mart=ensembl)
filters = listFilters(ensembl_h) # query # from
filters$name[grep('symbo', filters$name)] # hgnc_symbol

attributes = listAttributes(ensembl_h) # to
attributes$name[grep('name', attributes$name)] # wikigene_name external_gene_name
attributes$name[grep('hsa', attributes$name)] # 

geneDescr =  getBM(attributes=c('hgnc_symbol', 'description'), 
                   filters = 'hgnc_symbol', 
                   values = tolookagain$gene, 
                   mart = ensembl_h)
geneDescr2 = geneDescr %>%
  separate(description, into = c("description", "source"), sep = " \\[Source:", remove = F) %>%
  separate(source, into = c("source", "acc"), sep = ";Acc:", remove = T) %>%
  mutate(acc = gsub("\\]", "", acc),
         acc = gsub("HGNC:", "", acc))

go_geneDescr2 = tolookagain %>% #dplyr::select(tissue, gene, transcript) %>% 
  left_join(geneDescr2, by = c("gene" = "hgnc_symbol")) %>% 
  dplyr::select(-matches(".x")) %>% 
  dplyr::rename_with(., ~ (gsub(".y", "", .x)))
head(go_geneDescr2)

good2 = good %>% bind_rows(go_geneDescr2) 
dim(good2) # 310

path2 = paste0(output_path, "Hn_vs_Hg_TMM__Hn_list_biomart_", Sys.Date(), ".tsv")
# write_tsv(good2, path2, na= "")

good2_without_transcript = good2 %>% dplyr::select(-transcript) %>% distinct()
dim(good2_without_transcript) # 174
