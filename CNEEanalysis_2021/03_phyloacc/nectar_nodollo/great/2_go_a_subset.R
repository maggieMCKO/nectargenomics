# goal: for each KEGG term in our custom annotation, take genes, inner join the regulatory domain file by gene > output the subset BED file 

args = commandArgs(trailingOnly = TRUE)
nProc = args[1]
regDoms_path = args[2] # full path of great output
# regDoms_path = "/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/05_perm_phyloacc/1_rmPar/1_great/galGal6_gene_tss_great_basalPlusExtension.bed"

#### run GO analysis ####
# library(org.Gg.eg.db) # BiocManager::install("org.Gg.eg.db", lib)
library(tidyverse)
library(parallel)

set.seed(100)

project_path = "/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/03_phyloacc/2024ver/nectar_nodollo/great/"

# load custom annotation
path = paste0("/home/mpg08/mko/Nectar/analysis/CNEEanalysis_2021/00_inputs/
/validated.final.isoformes.galGal6.ncbi_appris_kegg_go_hsa_2022-11-21.tsv")
kegg_go_anno = read_tsv(path) 
go_annoC = kegg_go_anno %>% filter(!is.na(go_id)) %>%
  dplyr::select(gene, go_id, goTerm, ont) %>%
  filter(!is.na(gene)) %>%
  distinct() %>% 
  mutate(go_id = gsub(":", "_", go_id))
dim(go_annoC) # 325923      4
n_distinct(go_annoC$gene) # 15596
n_distinct(go_annoC$go_id) # 17917
# go_annoC %>% group_by(ont) %>% summarise(n = n_distinct(go_id)) # ignore NA
# 1 biological_process 11876
# 2 cellular_component  1783
# 3 molecular_function  4239
# 4 NA                    19
# x = go_annoC %>% filter(is.na(ont))
rm(kegg_go_anno)

go_anno_unique_terms = go_annoC %>% dplyr::select(go_id, ont) %>% distinct()

# load great regDom
regDoms = read_tsv(regDoms_path, col_names = F) %>% 
  dplyr::rename(gene = X4)
out_regDoms_basename = gsub(".bed", "", basename(regDoms_path))
out_regDoms_cleanname = gsub("galGal6_gene_tss_great_", "", out_regDoms_basename)

# create output dir
out_path = paste0(project_path, "/go/", out_regDoms_cleanname)
dir.create(out_path, recursive = TRUE)

# function to output bed file
func = function(term, in_ont){
  # term = go_anno_unique_terms_ont$go_id[1]
  print(paste0("current term: ", term))
  
  go_anno_sub = go_annoC %>% filter(go_id == term)
  regDoms_sub = regDoms %>% filter(gene %in% go_anno_sub$gene)
  
  if(nrow(regDoms_sub) >0){
    output_path = paste0(out_path, "/", in_ont, "/bed/", term, ".bed")
    write_tsv(regDoms_sub, file = output_path, na = "", col_names = FALSE)
  }
}

# run by ontogeny
onts = c("biological_process", "molecular_function", "cellular_component")
sapply(onts, function(ont_in){
  # ont_in = "molecular_function"
  print(paste0("current ont: ", ont_in))
  
  go_anno_unique_terms_ont = go_anno_unique_terms %>% filter(ont == ont_in) %>% 
    dplyr::select(-ont) %>% distinct()
  # output unique kegg pathways
  out_list_path = paste0(project_path, "/go/unique_go_", ont_in,".lst")
  write_tsv(go_anno_unique_terms_ont, file = out_list_path, col_names = FALSE)
  
  out_ont_path = paste0(project_path, "/go/", out_regDoms_cleanname, "/", ont_in, "/bed/")
  dir.create(out_ont_path, recursive = TRUE)
  
  term_list = go_anno_unique_terms_ont$go_id
  RunFunc = mclapply(X = term_list, FUN = func, mc.cores = nProc, in_ont = ont_in) 
  
})

