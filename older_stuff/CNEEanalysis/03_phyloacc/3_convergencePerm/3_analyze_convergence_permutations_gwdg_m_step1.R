#revised Oct 2018

library(tidyverse)
library(ape)

#functions
ctb <- function(x, cutoff = 0.90) {
  ifelse(x >= cutoff, 1, 0)
}

#load permutation results
wdir<-"/home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/3_convergencePerm"

runs<-list()
for (whichset in c("original")) {
  spec_patt<-glob2rx(paste0(whichset, "_perms_run*.tsv"))
  files<-list.files(path=paste0(wdir, "/convperms"), pattern=spec_patt, full.names = TRUE)
  results<-list()
  for (file in files) {
    results[[file]] <- read_tsv(file) 
  }
  runs[[whichset]] <- bind_rows(results, .id="file")
}

perms<-bind_rows(runs, .id="set") %>% select(set, version, test, tips, count) %>% distinct(set, version, test, tips, .keep_all = TRUE)

path = paste0(wdir, "/convperms/combined_", Sys.Date(), ".tsv")
write_tsv(perms, path)

