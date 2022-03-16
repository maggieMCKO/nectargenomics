library(tidyverse)


path = paste0( "../../galGal.gff.bed")
x = read_tsv(path, col_names = F) 
dim(x) # 17155 > rm chrM > 17142
x2 = x %>% mutate(gene = gsub("ID=gene-(\\w+-*\\.*\\w*);.*", "\\1", X10))
cc = sort(c(which(duplicated(x2$gene)), which(duplicated(x2$gene, fromLast = T))))
x2_dup = x2[cc,] # 0
length(unique(x2$gene)) == nrow(x) # TRUE => no duplicates
cc2 = sort(c(which(duplicated(x2$X4)), which(duplicated(x2$X4, fromLast = T))))
x2_dup = x2[cc2,] # 0
length(unique(x2$X4)) == nrow(x) # TRUE => no duplicates

path = paste0(project_path, "../../galGal.genes.bed")
pythonout = read_tsv(path, col_names = F) 
dim(pythonout) # 17155 > rm chrM > 17142
cc = sort(c(which(duplicated(pythonout$X4)), which(duplicated(pythonout$X4, fromLast = T))))
pythonout_dup = pythonout[cc,] # 28

x2_sub = lapply(1:nrow(pythonout_dup), function(i){
  b = pythonout_dup$X4[i]
  x2 %>% filter(grepl(b, X4)) %>% dplyr::select(all_of(c(1:4, 10, 11)))
}) %>% bind_rows() %>% distinct()
write_tsv(x2_sub, paste0(project_path, "../../duplicates.tsv"))