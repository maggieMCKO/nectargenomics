#revised Oct 2018

library(tidyverse)
library(ape)

# #functions
# ctb <- function(x, cutoff = 0.90) {
#   ifelse(x >= cutoff, 1, 0)
# }
# 
# #load permutation results
# project_path<-"/home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/3_convergencePerm"
# 
# runs<-list()
# for (whichset in c("original")) {
#   spec_patt<-glob2rx(paste0(whichset, "_perms_run*.tsv"))
#   files<-list.files(path=paste0(project_path, "/convperms"), pattern=spec_patt, full.names = TRUE)
#   results<-list()
#   for (file in files) {
#     results[[file]] <- read_tsv(file) 
#   }
#   runs[[whichset]] <- bind_rows(results, .id="file")
# }
# 
# perms<-bind_rows(runs, .id="set") %>% select(set, version, test, tips, count) %>% distinct(set, version, test, tips, .keep_all = TRUE)
# 
# path = paste0(project_path, "/convperms/combined_", Sys.Date(), ".tsv")
# # write_tsv(perms, path)


# in local PC
project_path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/3_convergencePerm/")
FontSize = 3
AxisTxFontSizeSize = 8
AxisTxFontSizeSize_s = 5
AxisTitleFontSizeSize = 10
Factor_mmtoin = 0.0393701
Width_HalfCol = 85*Factor_mmtoin 

openfile <- function(filepath){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  type = tolower(os);
  switch(os, Linux = system(paste0("xdg-open ", filepath)),
         Windows = system(paste0("open \"", filepath, "\"")),
         osx = system(paste0("open \"", filepath, "\""))
  )}

theme_m = 
  theme(
    plot.background = element_blank(),  #element_rect(fill = "transparent"), # defualt is white
    plot.title = element_text(hjust = 0.5, size = AxisTitleFontSizeSize),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(), #element_rect(color = "black", fill= NA, size = 0.1), # 
    panel.spacing.y = unit(0.2, "lines"),
    axis.line = element_line(color = "black", size = 0.1),
    axis.title.x = element_text(colour = "black", size = AxisTitleFontSizeSize),
    axis.title.y = element_text(colour = "black", size = AxisTitleFontSizeSize), # , hjust = 1, vjust = 0.5
    axis.text.x = element_text(colour = "black", size = AxisTxFontSizeSize),
    axis.text.y = element_text(colour = "black", size = AxisTxFontSizeSize), # , hjust = 1
    axis.ticks = element_line(colour = "black", size = 0.05),
    strip.text = element_text(colour = "black", size = AxisTxFontSizeSize, 
                              margin = margin(2,0,2,0, "pt")),
    legend.key = element_blank(),
    legend.text = element_text(colour = "black", size = AxisTxFontSizeSize),
    legend.title = element_text(colour = "black", size = AxisTitleFontSizeSize),
    legend.position = "none")

# path = paste0(project_path, "/convperms/combined_2020-07-21.tsv")
path = paste0(project_path, "/convperms/combined_2020-07-29.tsv")
perms = read_tsv(path)


#load real data
path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/1_runphyloacc/1_run_phyloacc/phyloacc_score_postZ.tsv")
# path = paste0("/home/mpg08/mko/Nectar/analysis/CNEEanalysis/03_phyloacc/1_run_phyloacc/phyloacc_score_postZ.tsv") # gwdg
cnee_orig_ori <- read_tsv(path)

cnee_orig = cnee_orig_ori %>% rename(cnee = ID ) %>% mutate(version = "basic") %>% 
  select(version, No., cnee, grep("_3", names(cnee_orig_ori))) # take acc _3 only
names(cnee_orig) = gsub("_3", "", names(cnee_orig))
rm(cnee_orig_ori)

path = paste0(project_path, "/input/nonconserved_4d_named_tree.nwk")
phy_orig <- read.tree(path)



#compute observed convergences with same approach as permutations

conv_real <- function(DF, targets, tips, cutoff = 0.90) {
  count_acc_targets <- function(tip, prob_acc, cutoff, targets) {
    targets_selected <- prob_acc[tip %in% targets]
    return(sum(targets_selected > cutoff))
  }
  
  #compute 
  # targets = target
  # DF = cnee_orig
  # tips = test_tips
  # cutoff = 0.9
  
  number<-length(targets)
  count<-DF %>% dplyr::select(cnee, any_of(tips)) %>% 
  # count<-DF %>% dplyr::select(cnee, one_of(tips)) %>% 
    gather(key = "tip", value = "prob_acc", -cnee) %>%
    group_by(cnee) %>%
    summarize(total_acc = sum(prob_acc > cutoff), target_acc = count_acc_targets(tip, prob_acc, cutoff, targets)) %>%
    dplyr::filter(target_acc == number & target_acc == total_acc) %>%
    tally() %>% pull(n)
  
  return(tibble(count=count, tips=paste0(targets, collapse="-")))
}

cross_real <- function(DF, targets, tips, cross_target, cutoff = 0.90) {
  count_acc_targets <- function(tip, prob_acc, cutoff, targets) {
    targets_selected <- prob_acc[tip %in% targets]
    return(sum(targets_selected > cutoff))
  }
  
  #compute 
  number<-length(targets)
  count<-DF %>% dplyr::select(cnee, one_of(tips)) %>% 
    gather(key = "tip", value = "prob_acc", -cnee) %>%
    group_by(cnee) %>%
    summarize(total_acc = sum(prob_acc > cutoff), 
              target_acc = count_acc_targets(tip, prob_acc, cutoff, targets), 
              cross_target_acc = count_acc_targets(tip, prob_acc, cutoff, cross_target)) %>%
    dplyr::filter(target_acc == number, target_acc == (total_acc-1), cross_target_acc == 1) %>%
    tally() %>% pull(n)
  
  return(tibble(count=count, tips=paste0(c(targets, cross_target), collapse="-")))
}



##
target = c("HLphyNov1", "HLlicCas1", "HLtriMol2", "HLcalAnn5", "HLfloFus1")
nontarget_tips_orig <- setdiff(phy_orig$tip.label, target)

test_tips <- phy_orig$tip.label

# select 3 species randomly ====
test_targets = combn(target, 3)
ratite_gm_dollo = lapply(1:ncol(test_targets), function(i){
  test_targets = test_targets[, i]
  test_tips <- c(test_targets, nontarget_tips_orig)
  cnee_orig %>% filter(version=="basic") %>%
    conv_real(test_targets, test_tips)
})

gm_dollo <- bind_rows(ratite_gm_dollo, .id="species")

p1 = perms %>% filter(version=="basic", set=="original", test=="conv_3") %>% 
  ggplot(aes(count)) + 
  geom_histogram(binwidth = 1) + 
  coord_cartesian(xlim=c(0,50)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_m + 
  
  # geom_segment(aes(x=mean(gm_dollo$count), xend=mean(gm_dollo$count), y=50, yend=0),
  #              arrow=arrow(length=unit(0.5, "cm")), colour="red") +
  geom_jitter(data=gm_dollo, aes(count, y=4), colour="red", width=0, height=3, size=2) +
  geom_text(data = subset(gm_dollo, su = count > 0), aes(label = tips, y = 50), size = FontSize, angle = 60, hjust = 0) +
  ylab("Count") + xlab("Number of Convergently Accelerated Elements"); p1


graph_path = paste0(project_path, "PermConvergence_", Sys.Date(), ".pdf")       # all data
pdf(graph_path, width = Width_HalfCol*1.5, height = Width_HalfCol*1, pointsize = AxisTxFontSizeSize, onefile = TRUE)
# png(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*0.625, unit = "in", res = 600)
p1
dev.off()
openfile(graph_path)

# select 3 clades randomly ====
honeyeaters = c("HLphyNov1", "HLlicCas1")
lori = "HLtriMol2"
hummingbirds = c("HLcalAnn5", "HLfloFus1")

test_targets2 = data.frame()
for (hum in hummingbirds) {
  for (honeyeater in honeyeaters){
    # test_targets = c(test_targets, list(lori, hum, honeyeater)) 
    x = data.frame(lori, hum, honeyeater)
    test_targets2 = test_targets2 %>% bind_rows(x)
  }
}
test_targets2 = t(test_targets2)

ratite_gm_dollo2 = lapply(1:ncol(test_targets2), function(i){
  test_targets = test_targets2[, i]
  test_tips <- c(test_targets, nontarget_tips_orig)
  cnee_orig %>% filter(version=="basic") %>%
    conv_real(test_targets, test_tips)
})

# ratite_gm_dollo2 = lapply(1:nrow(test_targets2), function(i){
#   test_targets = unlist(test_targets2[i, ] )
#   test_tips <- c(test_targets, nontarget_tips_orig)
#   cnee_orig %>% filter(version=="basic") %>%
#     conv_real(test_targets, test_tips)
# })

gm_dollo2 <- bind_rows(ratite_gm_dollo2, .id="species")

p2 = perms %>% filter(version=="basic", set=="original", test=="conv_3") %>% 
  ggplot(aes(count)) + 
  geom_histogram(binwidth = 1) + 
  coord_cartesian(xlim=c(0,50)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_m + 
  # geom_segment(aes(x=mean(gm_dollo$count), xend=mean(gm_dollo$count), y=50, yend=0),
  #              arrow=arrow(length=unit(0.5, "cm")), colour="red") +
  geom_jitter(data=gm_dollo2, aes(count, y=4), colour="red", width=0, height=3, size=2) +
  # geom_text(data = gm_dollo2, aes(label = tips, y = 50), size = FontSize, angle = 60, hjust = 0) +
  ylab("Count") + xlab("Number of Convergently Accelerated Elements"); p2


graph_path = paste0(project_path, "PermConvergence_clade_", Sys.Date(), ".pdf")       # all data
pdf(graph_path, width = Width_HalfCol*1.5, height = Width_HalfCol*1, pointsize = AxisTxFontSizeSize, onefile = TRUE)
# png(graph_path, width = Width_HalfCol*2, height = Width_HalfCol*0.625, unit = "in", res = 600)
p2
dev.off()
openfile(graph_path)

## tim's
# path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/2_GO/Tim/final_original_cnee.tsv.gz")
# cnee_orig_t <- read_tsv(path)

# path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/2_GO/Tim/cnees_extended.tsv.gz")
# cnee_ext_t <- read_tsv(path)

# path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis/03_phyloacc/3_convergencePerm/Tim/original.phy")
# phy_orig_t <- read.tree(path)te0(path_to_data, "/", "final_original_cnee.tsv.gz"))

# all_neo_orig <- cnee_orig_t %>% select(taeGut:anaPla) %>% names()
# all_neo_ext <- cnee_ext_t %>% select(taeGut:nanBra,uriPel:anaPla) %>% names()
# all_tin <- c("eudEle", "notPer", "tinGut", "cryCin")


#ratite dollo loop
ratite_gm_dollo <- list()
ratite_gl_dollo <- list()
for (third_sp in c("rheAme", "rhePen", "droNov", "casCas", "aptHaa", "aptOwe", "aptRow")) {
  test_targets <- c("strCam", "anoDid", third_sp)
  test_tips <- c(test_targets, all_neo_ext, all_tin)
  ratite_gm_dollo[[third_sp]] <- cnee_ext_t %>% filter(version=="gain") %>%
    conv_real(test_targets, test_tips)
  ratite_gl_dollo[[third_sp]] <- cnee_ext_t %>% filter(version=="gain_gap") %>%
    conv_real(test_targets, test_tips)
}

gm_dollo <- bind_rows(ratite_gm_dollo, .id="species")
gl_dollo <- bind_rows(ratite_gl_dollo, .id="species")

perms %>% filter(version=="gain", set=="original", test=="neo_conv_3") %>% ggplot(aes(count)) + 
  geom_histogram(binwidth = 1) + 
  coord_cartesian(xlim=c(0,50)) + theme_classic() + 
  geom_segment(aes(x=mean(gm_dollo$count), xend=mean(gm_dollo$count), y=50, yend=0), arrow=arrow(length=unit(0.5, "cm")), colour="red") +
  geom_jitter(data=gm_dollo, aes(count, y=4), colour="red", width=0, height=3, size=2) +
  ylab("Count") + xlab("Number of Convergently Accelerated Elements")

####

#alt versions
perms %>% filter(version=="gain_gap", set=="original", test=="neo_conv_3") %>% ggplot(aes(count)) + 
  geom_histogram(binwidth = 1) + 
  coord_cartesian(xlim=c(0,50)) + theme_classic() + 
  geom_segment(aes(x=mean(gl_dollo$count), xend=mean(gl_dollo$count), y=50, yend=0), arrow=arrow(length=unit(0.5, "cm")), colour="red") +
  geom_jitter(data=gl_dollo, aes(count, y=4), colour="red", width=0, height=3, size=2) +
  ylab("Count") + xlab("Number of Convergently Accelerated Elements")


#compute p-values
x=mean(gl_dollo$count)
perms %>% filter(version=="gain_gap", set=="original", test=="neo_conv_3") %>% 
  summarize(pval=(sum(count >= x)+1)/length(count))
x=median(gl_dollo$count)
perms %>% filter(version=="gain_gap", set=="original", test=="neo_conv_3") %>% 
  summarize(pval=(sum(count >= x)+1)/length(count))

x=mean(gm_dollo$count)
perms %>% filter(version=="gain", set=="original", test=="neo_conv_3") %>% 
  summarize(pval=(sum(count >= x)+1)/length(count))
x=median(gm_dollo$count)
perms %>% filter(version=="gain", set=="original", test=="neo_conv_3") %>% 
  summarize(pval=(sum(count >= x)+1)/length(count))

#numbers
mean(gl_dollo$count)
median(gl_dollo$count)
mean(gm_dollo$count)
median(gm_dollo$count)

#numbers from perms
perms %>% filter(test=="neo_conv_3") %>% with(., table(set, version))

perms %>% filter(test=="neo_conv_3", set=="original", version=="gain") %>% summarize(mean(count))

#counts / supplemental figures
#Fig S11 - plots of number of accelerated/convergent elements for gain, original dataset

cnee_gain_orig_conv <- cnee_orig %>% filter(dataset=="orig_v2_phyloAcc-gain", 
                                            logBF1 >= 10, logBF2 > 1,
                                            (it_pp_loss + ti_pp_loss) < 1, neo_tip_loss < 1) %>%
  mutate(floss_cl_bin = ctb(cd_pp_loss) + ctb(rh_pp_loss) + ctb(os_pp_loss) + ctb(ki_pp_loss) + ctb(mo_pp_loss),
         floss_cl_bin_dollo = ctb(ctb(cd_pp_loss) + ctb(rh_pp_loss) + ctb(ki_pp_loss)) + ctb(os_pp_loss) + ctb(mo_pp_loss)) %>%
  select(floss_cl_pp,floss_cl_pp_dollo,floss_cl_bin,floss_cl_bin_dollo,floss_sp_pp)

#different convergence metrics

#numbers

cnee_gain_orig_conv %>% summarize(count = sum(floss_cl_bin > 1.8))
cnee_gain_orig_conv %>% summarize(count = sum(floss_cl_bin_dollo > 1.8))
cnee_gain_orig_conv %>% summarize(count = sum(floss_cl_pp > 1.8))
cnee_gain_orig_conv %>% summarize(count = sum(floss_cl_pp_dollo > 1.8))

#percent
cnee_gain_orig_conv %>% summarize(count = sum(floss_cl_bin > 1.8)/length(floss_cl_bin))
cnee_gain_orig_conv %>% summarize(count = sum(floss_cl_bin_dollo > 1.8)/length(floss_cl_bin))
cnee_gain_orig_conv %>% summarize(count = sum(floss_cl_pp > 1.8)/length(floss_cl_bin))
cnee_gain_orig_conv %>% summarize(count = sum(floss_cl_pp_dollo > 1.8)/length(floss_cl_bin))


#plots

s10a <- cnee_gain_orig_conv %>% ggplot(aes(x=floss_cl_bin)) + geom_bar(fill="red") + theme_classic() + labs(x="Number of Accelerated Ratite Clades") + geom_vline(xintercept = 1.5)
s10b <- cnee_gain_orig_conv %>% ggplot(aes(x=floss_cl_bin_dollo)) + geom_bar(fill="red") + theme_classic() + labs(x="Number of Accelerated Ratite Clades") + geom_vline(xintercept = 1.5)
s10c<-cnee_gain_orig_conv %>% ggplot(aes(x=floss_cl_pp)) + geom_freqpoly(bins=50, col="red") + theme_classic() + labs(x="Posterior Estimate of Number of Independent Accelerations") + geom_vline(xintercept = 1.8)
s10d<-cnee_gain_orig_conv %>% ggplot(aes(x=floss_cl_pp_dollo)) + geom_freqpoly(bins=50, col="red") + theme_classic() + labs(x="Posterior Estimate of Number of Independent Accelerations") + geom_vline(xintercept = 1.8)
library(gridExtra)
grid.arrange(s10a,s10b,s10c,s10d,ncol=2)


#Fig S12 - consistency across models

cnee_orig %>% filter(logBF1 >= 10, logBF2 > 1, (it_pp_loss + ti_pp_loss) < 1) %>%
  ggplot(aes(floss_cl_pp, col=version)) + geom_freqpoly(bins=50, size=1.5) + theme_classic() +
  scale_color_brewer(palette = "Spectral")
ggsave("FigS12.pdf")


#Fig S13 - using gain version, compare original, extended, reduced

gain_comp_bf <- cnee_orig %>% filter(version=="gain") %>% select(cnee, logBF1, logBF2) %>% 
  inner_join(cnee_red %>% filter(version == "gain") %>% select(cnee, logBF1, logBF2), by=c("cnee" = "cnee"))

base_bf_plot <- gain_comp_bf %>% ggplot(aes(x=logBF1.x, y=logBF1.y)) + geom_hex(binwidth=c(3,3)) + theme_classic() + geom_hline(yintercept = 10, col="red") + geom_vline(xintercept = 10, col="red") + labs(x = "logBF1 (default dataset)", y="logBF1 (reduced dataset)") 

figS12A<-base_bf_plot + annotate("text", label = "Sig. in reduced, \n not in default: 256", x=-30, y=100) +
  annotate("text", label = "Sig. in both: 3527", x=200, y=100) +
  annotate("text", label = "Sig. in default, \n not in reduced: 1290", x=200, y=-25)

table(gain_comp_bf$logBF1.x >= 10, gain_comp_bf$logBF1.y >= 10)

gain_comp_pp <- cnee_orig %>% filter(version=="gain") %>% select(cnee, floss_cl_pp, logBF1, logBF2) %>% 
  full_join(cnee_red %>% filter(version == "gain") %>% select(cnee, floss_cl_pp, logBF1, logBF2), by=c("cnee" = "cnee")) %>% mutate(original = replace_na(floss_cl_pp.x, 0), reduced = replace_na(floss_cl_pp.y,0)) %>% 
  filter(logBF1.x >= 10 & logBF2.x >= 1 | logBF1.y >= 10 & logBF2.y >= 1)

figS12B<-gain_comp_pp %>% filter() %>% ggplot(aes(x=original, y=reduced)) + geom_point(alpha=0.2) + theme_classic()
grid.arrange(figS12A,figS12B)