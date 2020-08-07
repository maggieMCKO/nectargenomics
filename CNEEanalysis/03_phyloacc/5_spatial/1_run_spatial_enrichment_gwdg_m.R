#code updated for revisions Oct 2018
# [maggie] was analyze_spatial_enrichment_gwdg_m.R

# setwd("/n/holylfs/LABS/edwards_lab/tsackton/RATITE_PAPER_DATA_FREEZE/ratite-genomics/07_cnee_analysis/") #sorry!
library(tidyverse)
library(GenomicRanges) # bioconductor


# path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/02_phylop/galGal6_final_merged_CNEEs_named_fixchr_justchr.bed")
path = paste0("/home/mpg08/mko/Nectar/analysis/CNEEanalysis/02_phylop/galGal6_final_merged_CNEEs_named_fixchr_justchr.bed") # gwdg
pos_gg4_uscs<-read_tsv(path, col_names = c("chr", "start", "end", "cnee"))

chr_to_remove <- pos_gg4_uscs %>% count(chr) %>% filter(n < 150) %>% pull(chr)

# path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/1_runphyloacc/1_run_phyloacc/phyloacc_score_postZ.tsv")
path = paste0("/home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/1_run_phyloacc/phyloacc_score_postZ.tsv.gz") # gwdg
cnee_orig_ori <- read_tsv(path)
names(cnee_orig_ori) = gsub("ID", "cnee", names(cnee_orig_ori))
cnee_orig = cnee_orig_ori %>%
  dplyr::select(No., cnee, logBF1, logBF2,
                filter_non_target_5sp, filter_non_target_humm2anc, filter_non_target_humm2,
                filter_non_target_humm2anc_triMol, filter_non_target_humm2anc_honeyanc, filter_non_target_honey2anc_triMol) %>%
  full_join(pos_gg4_uscs, by=c("cnee" = "cnee")) %>%
  mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 , TRUE, FALSE)) %>%
  distinct(cnee, .keep_all=TRUE) %>%
  dplyr::select(chr, start, end, cnee, rar,
                filter_non_target_5sp, filter_non_target_humm2anc, filter_non_target_humm2,
                filter_non_target_humm2anc_triMol, filter_non_target_humm2anc_honeyanc, filter_non_target_honey2anc_triMol) %>%
  filter(!(chr %in% chr_to_remove)) %>%
  mutate(version = "basic")
rm(cnee_orig_ori)

# cnee_orig <- read_tsv("final_original_cnee.tsv.gz") %>%
#   select(version, cnee, logBF1, logBF2, it_pp_loss, ti_pp_loss, neo_tip_loss, floss_cl_pp, floss_cl_pp_dollo) %>%
#   full_join(pos_gg4_uscs, by=c("cnee" = "cnee")) %>%
#   mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 & (it_pp_loss + ti_pp_loss) < 1 & neo_tip_loss < 1, TRUE, FALSE),
#          crar = ifelse(rar & floss_cl_pp >= 1.8, TRUE, FALSE),
#          crar_dollo = ifelse(rar & floss_cl_pp_dollo >= 1.8, TRUE, FALSE)) %>%
#   distinct(cnee, version, .keep_all=TRUE) %>%
#   filter(!(chr %in% chr_to_remove)) %>%
#   select(chr, start, end, cnee, version, rar, crar, crar_dollo)


#note in ext2, convergence defined as ratites + cormorants
#cnee_ext2 <- read_tsv("final_extended_cnee.tsv.gz") %>%
#  select(version, cnee, logBF1, logBF2, it_pp_loss, ti_pp_loss, neo_tip_loss, #floss_cl_pp, floss_cl_pp_dollo, gc_pp_loss) %>%
#  full_join(pos_gg4_uscs, by=c("cnee" = "cnee")) %>%
#  mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 & (it_pp_loss + ti_pp_loss) < 1 & neo_tip_loss < 1, TRUE, FALSE),
#         crar = ifelse(rar & floss_cl_pp >= 1.8 & gc_pp_loss > 0.90, TRUE, FALSE),
#         crar_dollo = ifelse(rar & floss_cl_pp_dollo >= 1.8 & gc_pp_loss > 0.90, #TRUE, FALSE)) %>%
#  distinct(cnee, version, .keep_all=TRUE) %>%
#  filter(!(chr %in% chr_to_remove)) %>%
#  select(chr, start, end, cnee, version, rar, crar, crar_dollo)


#we'll use pbinom to compute a simple test of whether, given random distribution of RARs/cRARS, we'd expect X or more in a window
#q = number RARs in window
#size = number CNEEs in window
#prob = total RARs / total CNEEs (in set)
#lower.tail = F will give P-value of q or greater observed given prob

#need to define windows, then get q and size for each window, and compute p-value

#read in seqlengths
# path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/0_preprocessing/1_split_sp/galGal6chromSize.txt")
path = paste0("/home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/0_preprocessing/1_way2/galGal6chromSize.txt") # gwdg
seqlen_df<-read_tsv(path, col_names = F) # %>% filter(X1 %in% chr_to_remove ==FALSE)
seqlen<-seqlen_df$X2
names(seqlen)<-seqlen_df$X1

#function to do analysis

compute_spatial_results <- function(DF, outname) {
  
  # DF = cnee_orig
  # dir.create(dirname(outname))
  
  spa_res <- list()
  
  for (ver in c("basic")) {
    
    cnee <- DF %>% filter(version == ver)
    
    #convert to window
    cnee_ranges<-makeGRangesFromDataFrame(cnee, keep.extra.columns = TRUE, ignore.strand=TRUE, seqinfo = seqlen)
    windows_500kb<-tileGenome(seqinfo(cnee_ranges), tilewidth=500000, cut.last.tile.in.chrom = TRUE)
    windows_1000kb<-tileGenome(seqinfo(cnee_ranges), tilewidth=1000000, cut.last.tile.in.chrom = TRUE)
    windows_100kb<-tileGenome(seqinfo(cnee_ranges), tilewidth=100000, cut.last.tile.in.chrom = TRUE)
    
    #1Mb sliding window
    #hacky
    windows_chr<-tileGenome(seqinfo(cnee_ranges), tilewidth=max(seqlengths(cnee_ranges)), cut.last.tile.in.chrom=TRUE)
    windows_1Mb_slide <- slidingWindows(windows_chr, 1e6, 100000) %>% unlist
    
    rar_prob<-sum(cnee_ranges$rar, na.rm=T)/length(cnee_ranges$rar)
    filter_non_target_5sp_prob<-sum(cnee_ranges$filter_non_target_5sp, na.rm=T)/length(cnee_ranges$filter_non_target_5sp)
    filter_non_target_humm2anc_prob<-sum(cnee_ranges$filter_non_target_humm2anc, na.rm=T)/length(cnee_ranges$filter_non_target_humm2anc)
    filter_non_target_humm2_prob<-sum(cnee_ranges$filter_non_target_humm2, na.rm=T)/length(cnee_ranges$filter_non_target_humm2)
    filter_non_target_humm2anc_triMol_prob<-sum(cnee_ranges$filter_non_target_humm2anc_triMol, na.rm=T)/length(cnee_ranges$filter_non_target_humm2anc_triMol)
    filter_non_target_humm2anc_honeyanc_prob<-sum(cnee_ranges$filter_non_target_humm2anc_honeyanc, na.rm=T)/length(cnee_ranges$filter_non_target_humm2anc_honeyanc)
    filter_non_target_honey2anc_triMol_prob<-sum(cnee_ranges$filter_non_target_honey2anc_triMol, na.rm=T)/length(cnee_ranges$filter_non_target_honey2anc_triMol)
    # crar_dollo_prob<-sum(cnee_ranges$crar_dollo, na.rm=T)/length(cnee_ranges$crar_dollo)

    get_stats_for_window <- function(query_ranges, window_to_test) {
      window_df <- as.data.frame(window_to_test)
      stopifnot(length(window_df$seqnames)==1)
      subsetByOverlaps(query_ranges, window_to_test) %>%
        as.data.frame %>%
        summarize(rar_n = sum(rar, na.rm=T),
                  filter_non_target_5sp_n = sum(filter_non_target_5sp, na.rm=T),
                  filter_non_target_humm2anc_n = sum(filter_non_target_humm2anc, na.rm=T),
                  filter_non_target_humm2_n = sum(filter_non_target_humm2, na.rm=T),
                  filter_non_target_humm2anc_triMol_n = sum(filter_non_target_humm2anc_triMol, na.rm=T),
                  filter_non_target_humm2anc_honeyanc_n = sum(filter_non_target_humm2anc_honeyanc, na.rm=T),
                  filter_non_target_honey2anc_triMol_n = sum(filter_non_target_honey2anc_triMol, na.rm=T),
                  # crar_dollo_n = sum(crar_dollo, na.rm=T),
                  total_n=n()) %>%
        mutate(rar_pval = pbinom(q=rar_n, size=total_n, prob=rar_prob, lower.tail=FALSE) +
                 dbinom(x=rar_n, size=total_n, prob=rar_prob),
               filter_non_target_5sp_pval = pbinom(q=filter_non_target_5sp_n, size=total_n, prob=filter_non_target_5sp_prob, lower.tail=FALSE) +
                 dbinom(x=filter_non_target_5sp_n, size=total_n, prob=filter_non_target_5sp_prob),
               filter_non_target_humm2anc_pval = pbinom(q=filter_non_target_humm2anc_n, size=total_n, prob=filter_non_target_humm2anc_prob, lower.tail=FALSE) +
                 dbinom(x=filter_non_target_humm2anc_n, size=total_n, prob=filter_non_target_humm2anc_prob),
               filter_non_target_humm2_pval = pbinom(q=filter_non_target_humm2_n, size=total_n, prob=filter_non_target_humm2_prob, lower.tail=FALSE) +
                 dbinom(x=filter_non_target_humm2_n, size=total_n, prob=filter_non_target_humm2_prob),
               filter_non_target_humm2anc_triMol_pval = pbinom(q=filter_non_target_humm2anc_triMol_n,
                                                               size=total_n, prob=filter_non_target_humm2anc_triMol_prob, lower.tail=FALSE) +
                 dbinom(x=filter_non_target_humm2anc_triMol_n, size=total_n, prob=filter_non_target_humm2anc_triMol_prob),
               filter_non_target_humm2anc_honeyanc_pval = pbinom(q=filter_non_target_humm2anc_honeyanc_n,
                                                                 size=total_n, prob=filter_non_target_humm2anc_honeyanc_prob, lower.tail=FALSE) +
                 dbinom(x=filter_non_target_humm2anc_honeyanc_n, size=total_n, prob=filter_non_target_humm2anc_honeyanc_prob),
               filter_non_target_honey2anc_triMol_pval = pbinom(q=filter_non_target_honey2anc_triMol_n,
                                                                size=total_n, prob=filter_non_target_honey2anc_triMol_prob, lower.tail=FALSE) +
                 dbinom(x=filter_non_target_honey2anc_triMol_n, size=total_n, prob=filter_non_target_honey2anc_triMol_prob),
               # crar_dollo_pval = pbinom(q=crar_dollo_n, size=total_n, prob=crar_dollo_prob, lower.tail=FALSE) +
                 # dbinom(x=crar_dollo_n, size=total_n, prob=crar_dollo_prob),
               chr= window_df$seqnames, start=window_df$start, end=window_df$end) %>%
        select(chr, start, end, rar_n,
               filter_non_target_5sp_n, filter_non_target_humm2anc_n, filter_non_target_humm2_n,
               filter_non_target_humm2anc_triMol_n, filter_non_target_humm2anc_honeyanc_n, filter_non_target_honey2anc_triMol_n,
               # crar_n, crar_dollo_n,
               total_n, rar_pval,
               filter_non_target_5sp_pval, filter_non_target_humm2anc_pval, filter_non_target_humm2_pval,
               filter_non_target_humm2anc_triMol_pval, filter_non_target_humm2anc_honeyanc_pval, filter_non_target_honey2anc_triMol_pval)
               # crar_pval, crar_dollo_pval)
    }
    
    se_500kb<-lapply(1:length(windows_500kb), function(x) get_stats_for_window(cnee_ranges, windows_500kb[x])) %>% dplyr::bind_rows(.id="window") %>%
      filter(total_n > 0)
    se_500kb$rar_qval <- p.adjust(se_500kb$rar_pval, "fdr")
    se_500kb$filter_non_target_5sp_qval <- p.adjust(se_500kb$filter_non_target_5sp_pval, "fdr")
    se_500kb$filter_non_target_humm2anc_qval <- p.adjust(se_500kb$filter_non_target_humm2anc_pval, "fdr")
    se_500kb$filter_non_target_humm2_qval <- p.adjust(se_500kb$filter_non_target_humm2_pval, "fdr")
    se_500kb$filter_non_target_humm2anc_triMol_qval <- p.adjust(se_500kb$filter_non_target_humm2anc_triMol_pval, "fdr")
    se_500kb$filter_non_target_humm2anc_honeyanc_qval <- p.adjust(se_500kb$filter_non_target_humm2anc_honeyanc_pval, "fdr")
    se_500kb$filter_non_target_honey2anc_triMol_qval <- p.adjust(se_500kb$filter_non_target_honey2anc_triMol_pval, "fdr")
    # se_500kb$crar_qval <- p.adjust(se_500kb$crar_pval, "fdr")
    # se_500kb$crar_dollo_qval <- p.adjust(se_500kb$crar_dollo_pval, "fdr")
    se_500kb$window_size = "500kb"
    
    se_100kb<-lapply(1:length(windows_100kb), function(x) get_stats_for_window(cnee_ranges, windows_100kb[x])) %>% dplyr::bind_rows(.id="window") %>%
      filter(total_n > 0)
    se_100kb$rar_qval <- p.adjust(se_100kb$rar_pval, "fdr")
    se_100kb$filter_non_target_5sp_qval <- p.adjust(se_100kb$filter_non_target_5sp_pval, "fdr")
    se_100kb$filter_non_target_humm2anc_qval <- p.adjust(se_100kb$filter_non_target_humm2anc_pval, "fdr")
    se_100kb$filter_non_target_humm2_qval <- p.adjust(se_100kb$filter_non_target_humm2_pval, "fdr")
    se_100kb$filter_non_target_humm2anc_triMol_qval <- p.adjust(se_100kb$filter_non_target_humm2anc_triMol_pval, "fdr")
    se_100kb$filter_non_target_humm2anc_honeyanc_qval <- p.adjust(se_100kb$filter_non_target_humm2anc_honeyanc_pval, "fdr")
    se_100kb$filter_non_target_honey2anc_triMol_qval <- p.adjust(se_100kb$filter_non_target_honey2anc_triMol_pval, "fdr")
    # se_100kb$crar_qval <- p.adjust(se_100kb$crar_pval, "fdr")
    # se_100kb$crar_dollo_qval <- p.adjust(se_100kb$crar_dollo_pval, "fdr")
    se_100kb$window_size = "100kb"
    
    
    se_1000kb<-lapply(1:length(windows_1000kb), function(x) get_stats_for_window(cnee_ranges, windows_1000kb[x])) %>% dplyr::bind_rows(.id="window") %>%
      filter(total_n > 0)
    se_1000kb$rar_qval <- p.adjust(se_1000kb$rar_pval, "fdr")
    se_1000kb$filter_non_target_5sp_qval <- p.adjust(se_1000kb$filter_non_target_5sp_pval, "fdr")
    se_1000kb$filter_non_target_humm2anc_qval <- p.adjust(se_1000kb$filter_non_target_humm2anc_pval, "fdr")
    se_1000kb$filter_non_target_humm2_qval <- p.adjust(se_1000kb$filter_non_target_humm2_pval, "fdr")
    se_1000kb$filter_non_target_humm2anc_triMol_qval <- p.adjust(se_1000kb$filter_non_target_humm2anc_triMol_pval, "fdr")
    se_1000kb$filter_non_target_humm2anc_honeyanc_qval <- p.adjust(se_1000kb$filter_non_target_humm2anc_honeyanc_pval, "fdr")
    se_1000kb$filter_non_target_honey2anc_triMol_qval <- p.adjust(se_1000kb$filter_non_target_honey2anc_triMol_pval, "fdr")
    # se_1000kb$crar_qval <- p.adjust(se_1000kb$crar_pval, "fdr")
    # se_1000kb$crar_dollo_qval <- p.adjust(se_1000kb$crar_dollo_pval, "fdr")
    se_1000kb$window_size = "1000kb"
    
    #sliding windows
    se_slide<-lapply(1:length(windows_1Mb_slide), function(x) get_stats_for_window(cnee_ranges, windows_1Mb_slide[x])) %>%
      dplyr::bind_rows(.id="window") %>%
      filter(total_n > 0)
    se_slide$rar_qval <- p.adjust(se_slide$rar_pval, "fdr")
    se_slide$filter_non_target_5sp_qval <- p.adjust(se_slide$filter_non_target_5sp_pval, "fdr")
    se_slide$filter_non_target_humm2anc_qval <- p.adjust(se_slide$filter_non_target_humm2anc_pval, "fdr")
    se_slide$filter_non_target_humm2_qval <- p.adjust(se_slide$filter_non_target_humm2_pval, "fdr")
    se_slide$filter_non_target_humm2anc_triMol_qval <- p.adjust(se_slide$filter_non_target_humm2anc_triMol_pval, "fdr")
    se_slide$filter_non_target_humm2anc_honeyanc_qval <- p.adjust(se_slide$filter_non_target_humm2anc_honeyanc_pval, "fdr")
    se_slide$filter_non_target_honey2anc_triMol_qval <- p.adjust(se_slide$filter_non_target_honey2anc_triMol_pval, "fdr")
    # se_slide$crar_qval <- p.adjust(se_slide$crar_pval, "fdr")
    # se_slide$crar_dollo_qval <- p.adjust(se_slide$crar_dollo_pval, "fdr")
    se_slide$window_size = "1000kb_100kb_slide"
    
    spa_res[[ver]] <- rbind(se_100kb, se_500kb, se_1000kb, se_slide) %>% as.tibble
  }

  bind_rows(spa_res, .id="version") %>% write_tsv(outname)
}

outname_input = paste0("/home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/5_spatial/original_spatial_results.tsv")
compute_spatial_results(cnee_orig, outname_input)
# compute_spatial_results(cnee_ext, "extended_spatial_results.tsv")
# compute_spatial_results(cnee_red, "reduced_spatial_results.tsv")

