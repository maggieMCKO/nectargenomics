#revised Oct 2018 for paper revisions
library(tidyverse)

#loop to compute ecdfs for each run (extended, original, reduced, extended cormorant and version galgal4, galgal5)

wdir<-"/home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/2_GO/1_enrichment"

for (whichset in c("original")) {
  for (whichgenome in c("galgal6")) {
    for (whichont in c("BP", "MF")) {
    	spec_patt<-glob2rx(paste0(whichset, "_GO_", whichgenome, "_run*_", whichont, "_perm.tsv"))
    	files<-list.files(path=paste0(wdir, "/goperms"), pattern=spec_patt, full.names = TRUE)
    	results<-list()
    	for (file in files) {
    	  results[[file]] <- read_tsv(file) %>% filter(!is.na(target_frac) & target_frac < 1) %>% mutate(enrich = log2(target_frac/bg_frac))
    	}
    	tmp = bind_rows(results)
    	write_csv(tmp, paste0(wdir, "/goperms/", whichset, "_", whichgenome, "_", whichont, "_raw_", Sys.Date(), ".csv"))
    	
    	tmp %>% 
    	  group_by(version, set, ID) %>% 
    	  summarize(ecdf_frac = list(ecdf(target_frac)),
    	            ecdf_enrich = list(ecdf(enrich)),
    	            ecdf_logp = list(ecdf(logp.perm))) %>% 
    	  saveRDS(file=paste0(wdir, "/goperms/", whichset, "_", whichgenome, "_", whichont, "_", Sys.Date(), ".robj"))
    	}
  }
}
