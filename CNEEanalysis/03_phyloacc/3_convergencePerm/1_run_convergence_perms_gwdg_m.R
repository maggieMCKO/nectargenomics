#this code runs the permutations to test whether convergence is greater than expected
library(tidyverse)
library(parallel)
library(ape)

#setwd("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/")

args <- commandArgs(trailingOnly = TRUE)

# OUTNUM <- args[1] # never used
CORES <- args[2]
NPERM <- args[3]
path_to_data <-args[4] # == project_path

# project_path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/3_convergencePerm/")
# project_path = paste0("/home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/3_convergencePerm/") # gwdg
project_path = path_to_data

# test
# OUTNUM <- paste0(project_path, "/output/run", 1)
# CORES <- 1
# NPERM <- 5
# path_to_data <-project_path

## LOAD DATA ##

# path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/1_runphyloacc/1_run_phyloacc/phyloacc_score_postZ.tsv")
path = paste0("/home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/1_run_phyloacc/phyloacc_score_postZ.tsv") # gwdg
cnee_orig_ori <- read_tsv(path)
cnee_orig = cnee_orig_ori %>% rename(cnee = ID ) %>% mutate(version = "basic") %>% 
  select(version, No., cnee, grep("_3", names(cnee_orig_ori))) # take acc only
names(cnee_orig) = gsub("_3", "", names(cnee_orig))

# path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/3_convergencePerm/input/nonconserved_4d_named_tree.nwk")
path = paste0(project_path, "/input/nonconserved_4d_named_tree.nwk") # gwdg
phy_orig <- read.tree(path)
# phy_orig <- read.tree(paste0(path_to_data, "/", "original.phy"))


## tim's
# path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/2_GO/Tim/final_original_cnee.tsv.gz")
# cnee_orig_t <- read_tsv(path)
# # cnee_orig <- read_tsv(paste0(path_to_data, "/", "final_original_cnee.tsv.gz"))
# 
# path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/3_convergencePerm/Tim/original.phy")
# phy_orig_t <- read.tree(path)


## FUNCTIONS ##

#functions
tto <- function(x) {
  ifelse(x > 1, 1, x)
}

ctb <- function(x, cutoff=0.90, lower=FALSE) {
  if (lower) {
    ifelse(x <= cutoff, 1, 0)
  }
  else {
    ifelse(x >= cutoff, 1, 0)
  }
}

get_max <- function(x, cutoff = 0.90) {
  ifelse(max(x) > cutoff, 1, 0)
}

count_sister_taxa <- function(tree, tips) {
  #set up pairs
  pairs<-combn(tips, 2, simplify = FALSE)
  lapply(pairs, is.monophyletic, phy=tree) %>% unlist %>% sum
}

conv_sample <- function(perm=1, DF, number, tips, phy, cutoff = 0.90) {
  # DF = cnee_orig
  # number = 5
  # tips = nontarget_tips_orig
  # phy = phy_orig
  # cutoff = 0.9
  
  # Tim's
  # DF = cnee_orig_t
  # number = 3
  # phy = phy_orig_t
  # tips = neo_tips_orig
  # cutoff = 0.9
  
  num_sister = 10
  while (num_sister > 0) {
    targets<-sort(sample(tips, number))
    num_sister = count_sister_taxa(phy, targets); # print(num_sister)
  }
  #targets now has N (= number) random non-sister tips
  
  count_acc_targets <- function(tip, prob_acc, cutoff, targets) {
    targets_selected <- prob_acc[tip %in% targets]
    return(sum(targets_selected > cutoff))
  }
  
  #compute 
  
  # count<- DF %>% dplyr::select(cnee, one_of(tips)) %>% 
  count<- DF %>% dplyr::select(cnee, any_of(tips)) %>% 
    gather(key = "tip", value = "prob_acc", -cnee) %>%
    group_by(cnee) %>%
    summarize(total_acc = sum(prob_acc > cutoff), target_acc = count_acc_targets(tip, prob_acc, cutoff, targets)) %>%
    dplyr::filter(target_acc == number & target_acc == total_acc) %>%
    tally() %>% pull(n)
  
  return(tibble(count=count, tips=paste0(targets, collapse="-")))
} # conv_sample

conv_sample_cross <- function(perm=1, DF, number, tips, cross_tips, phy, cutoff = 0.90) {
  #bad progamming practice copying functions from above...oh well
  if (number > 1) {
    num_sister = 10
    while (num_sister > 0) {
      targets<-sort(sample(tips, number))
      num_sister = count_sister_taxa(phy, targets)  
    } 
  } else {
    targets<-sort(sample(tips, 1))
  }
  #targets now has N (= number) random non-sister tips
  
  count_acc_targets <- function(tip, prob_acc, cutoff, targets) {
    targets_selected <- prob_acc[tip %in% targets]
    return(sum(targets_selected > cutoff))
  }
  
  #cross target is the single "extra" species to compare
  cross_target<-sample(cross_tips, 1)
  
  #compute 
  
  count<-DF %>% dplyr::select(cnee, one_of(tips), cross_target) %>% 
    gather(key = "tip", value = "prob_acc", -cnee) %>%
    group_by(cnee) %>%
    summarize(total_acc = sum(prob_acc > cutoff), 
              target_acc = count_acc_targets(tip, prob_acc, cutoff, targets), 
              cross_target_acc = count_acc_targets(tip, prob_acc, cutoff, cross_target)) %>%
    dplyr::filter(target_acc == number, target_acc == (total_acc-1), cross_target_acc == 1) %>%
    tally() %>% pull(n)
  
  return(tibble(count=count, tips=paste0(c(targets, cross_target), collapse="-")))  
} # conv_sample_cross


#EACH DATASET NEEDS DIFFERENT PROCESSING SO NEED TO DO IN SEQUENCE

##ORIGINAL (+moa, -corm)

orig_perms<-list()
# red_perms<-list()
# ext_perms<-list()

target = c("HLphyNov1", "HLlicCas1", "HLtriMol2", "HLcalAnn5", "HLfloFus1")
nontarget_tips_orig <- setdiff(phy_orig$tip.label, target)
# neo_tips_orig <- phy_orig_t$tip.label[1:23]  # everything non-target
# neo_tips_ext <- phy_ext$tip.label[c(1:14,16:27)]
# tin_tips <- phy_ext$tip.label[35:38]

orig_perms[[1]] <- cnee_orig %>% filter(version=="basic") %>% 
  mclapply(1:NPERM, conv_sample, DF=., number=5, tips=nontarget_tips_orig, phy=phy_orig, mc.preschedule = TRUE, mc.cores = CORES) %>%
  bind_rows(.id="perm") %>%
  mutate(version = "basic", dataset="original", test="conv_5") %>% 
  distinct(tips, .keep_all=TRUE)

orig_perms[[2]] <- cnee_orig %>% filter(version=="basic") %>% 
  mclapply(1:NPERM, conv_sample, DF=., number=3, tips=nontarget_tips_orig, phy=phy_orig, mc.preschedule = TRUE, mc.cores = CORES) %>%
  bind_rows(.id="perm") %>%
  mutate(version = "basic", dataset="original", test="conv_3") %>% 
  distinct(tips, .keep_all=TRUE)

out_dir = paste0(path_to_data, "/", "convperms/")
dir.create(out_dir)

outpath = paste0(out_dir, "/original_perms_run", args[1], ".tsv")
orig_perms %>% bind_rows() %>% write_tsv(outpath)
