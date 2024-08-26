# goal: for each KEGG term in our custom annotation, take genes, inner join the regulatory domain file by gene > output the subset BED file 

args = commandArgs(trailingOnly = TRUE)
input1path = args[1] # ori(file to be subset) path
input2path = args[2] # file with biotype
outputpath = args[3] # output path

# input1path = "/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/05_perm_phyloacc/1_rmPar/1_great/galGal6_gene_tss.bed"
# input2path = "/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/05_perm_phyloacc/1_rmPar/1_great/galGal.genes.bed"
# outputpath = "/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/05_perm_phyloacc/1_rmPar/1_great/galGal6_gene_tss_proteincoding.bed"


#### run GO analysis ####
# library(org.Gg.eg.db) # BiocManager::install("org.Gg.eg.db", lib)
library(tidyverse)

input1 = read_tsv(input1path, col_names = F)
dim(input1) # 24108     4

input2 = read_tsv(input2path, col_names = F)
dim(input2) # way2: 24108     4; way3:  17485     9
input2_cleaned = input2 %>% 
  separate(X9, into = c("gene", "rest"), sep = ";", remove = F) %>% 
  mutate(gene = gsub("ID=gene-", "", gene))  %>% 
  separate(X9, into = c("rest", "biotype"), sep = "gene_biotype=", remove = F)  %>% 
  separate(biotype, into = c("biotype2", "rest2"), sep = ";", remove = F) %>% 
  dplyr::select(X1, X3, gene, biotype2) %>% 
  dplyr::rename(biotype = biotype2)
dim(input2_cleaned) #  17485     4
unique(input2_cleaned$biotype) # "protein_coding"

# gplots::venn(list(input1$X4, input2_cleaned$gene))
# onlya:intersect:onlyb = 6636:17472:13
# onlyinput1 = input1 %>% filter(!X4 %in% input2_cleaned$gene) 
# onlyx_cleaned2 = input2_cleaned %>% filter(!gene %in% input1$X4) # chrMT


# input2_prot = input2 %>% filter(X8 == "protein_coding")
# n_distinct(input2_prot$X6) # 11603

input1_sub = input1 %>% filter(X4 %in% input2_cleaned$gene)
n_distinct(input1_sub$X4) # wah1:11603>way3>17472

write_tsv(input1_sub, file = outputpath, col_names = F)

# check and some scrabbles below
# input1_notcoding = input1 %>% filter(!X4 %in% input2_cleaned$gene)

# input1_notcoding = input1 %>% filter(!X4 %in% input2_prot$X6)

# gplots::venn(list(input1$X4, input2$X6))
# all input2 are in input1, 6598 only in input1 when using way1 to generate input2, way2 output total over
# input1_only = input1 %>% filter(!X4 %in% input2$X6)

# gffpath = "/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/00_inputs/GCF_000002315.6_GRCg6a_genomic_chrtabFixed_justchr.gff"
# gff = read_tsv(gffpath, col_names = F)
# dim(gff) # 1589911      13
# gff_proteincoding = gff %>% filter(grepl('protein_coding', X9)) %>%
#   separate(X9, into = c("gene", "rest"), sep = ";", remove = F) %>% 
#   mutate(gene = gsub("ID=gene-", "", gene))
