library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
#args 1 is where the dir with pcoc output

project_path = args[1]
# project_path = "/home/mpg08/mko/Nectar/analysis/rnaseq/00_inputs/Kallisto_quntification_all_lib" 

target_files = list.files(project_path, "gene_abundance.tsv", full.names = TRUE, recursive = TRUE)

out = lapply(target_files, function(s){
  # s = '/home/mpg08/mko/Nectar/analysis/pcoc/00_inputs/pcoc_sim/s2/1500/RUN_20220516_114842/Tree_1/BenchmarkResults.tsv'
  
  tmp_path = s
  tmp_dir = dirname(tmp_path)
  tmp = unlist(strsplit(gsub(project_path, "", tmp_dir), "/"))
  # print(tmp)
  tmp_sample_id = tmp[2]
  tb = read_tsv(tmp_path, col_names = c("gene", "counts")) %>% bind_cols("sample_id" = tmp_sample_id)
  
}) %>% bind_rows()

path = paste0(project_path, "/../gene_abundance_all.tsv")
write_tsv(out, path, na = "")