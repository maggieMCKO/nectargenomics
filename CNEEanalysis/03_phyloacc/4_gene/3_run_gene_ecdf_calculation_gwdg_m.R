#revised Oct 2018 for paper revisions
library(tidyverse)

#loop to compute ecdfs for each run (extended, original, reduced, extended cormorant and version galgal4, galgal5)

wdir<-"/home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/4_gene"

for (whichset in c("original")) {
  for (whichgenome in c("galgal6")) {
    spec_patt<-glob2rx(paste0(whichset, "_gene_", whichgenome, "_run*_perm.tsv"))
    files<-list.files(path=paste0(wdir, "/geneperms"), pattern=spec_patt, full.names = TRUE)
    results<-list()
    for (file in files) {
      results[[file]] <- read_tsv(file) 
    }
    bind_rows(results) %>% group_by(version, set, gene) %>% 
      summarize(ecdf_gene = list(ecdf(rand_TRUE))) %>% 
      saveRDS(file=paste0(wdir, "/geneperms/", whichset, "_", whichgenome, ".robj"))
  }
}
