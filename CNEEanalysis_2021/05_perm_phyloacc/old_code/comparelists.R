library(tidyverse)
library(UpSetR)
library(biomaRt) # bioconductor

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

FontSize = 2.5
AxisTxFontSizeSize = 8
AxisTxFontSizeSize_s = 6
AxisTitleFontSizeSize = 10
Factor_mmtoin = 0.0393701
Width_HalfCol = 85*Factor_mmtoin 


project_path = paste0(getwd(), "/Seq_Data/globus/CNEEanalysis_2021/05_perm_phyloacc/1_prep_new/postPhyloAcc/")

absrel = read_tsv(paste0(getwd(), "/../../Presentation/Lab meeting/important_sldies/Results_FrKatya/20220321_absrelgenelists/under_selection_per_clade_0.05.noReg.tsv"))
unique(absrel$group)
absrel = absrel %>% 
  mutate(group = gsub("hmmbrds", "hummingbirds", group),
         group = gsub("honeyeaters", "honeyeaters_pardalote", group),
         group = gsub("nectar_parrots", "parrots", group),
         group = gsub("sunbirds", "sunbirds_flowerpecker", group)) %>%
  bind_cols("analysis" = "absrel")

spatial = read_tsv(paste0(project_path, "spatial/AccCneeEnrichedWindows_proteincoding_padj0.05.tsv")) %>%
  bind_cols("analysis" = "spatial")
unique(spatial$group)

gene = read_tsv(paste0(project_path, "gene/AccCneeEnrichedProteinCodingGenes_padj0.05.tsv")) %>%
  bind_cols("analysis" = "gene")
unique(spatial$group)



comb = absrel %>% bind_rows(spatial) %>% bind_rows(gene) %>%
  mutate(set = paste0(group, "_", analysis)) %>%
  dplyr::select(gene, set) %>%
  distinct() %>%
  mutate(val = 1)
length(unique(comb$set)) # 10
table(comb$set)
nrow(comb) # 1561
length(unique(comb$gene)) # 1013

# upsetR ####

## to wide
present = function(s){ifelse(!is.na(s), 1, 0)}
# x = sapply(comb$val, function(ss){present(ss)})
comb_wide = comb %>% 
  pivot_wider(id_cols = gene, names_from = set, values_from = val, values_fn = present, values_fill = 0) %>%
  as.data.frame()

unique(comb$set)

p = upset(comb_wide, 
          sets = rev(c(
            "hummingbirds_absrel",
            "hummingbirds_gene",
            "hummingbirds_spatial",
            
            "parrots_absrel",
            "honeyeaters_pardalote_absrel",
            "honeyeaters_pardalote_gene",
            "honeyeaters_pardalote_spatial",
            
            "sunbirds_flowerpecker_absrel",
            "sunbirds_flowerpecker_gene",
            "sunbirds_flowerpecker_spatial"
            )), keep.order = TRUE,
          main.bar.color = "#56B4E9",
          matrix.color = "#56B4E9",
          sets.bar.color = "#56B4E9",
          mb.ratio = c(0.6, 0.4),
          order.by = c("degree", "freq")); p


graph_path = paste0(project_path, "upset_absrel_spatial_gene_", Sys.Date(),".pdf"); graph_path
pdf(graph_path, width = Width_HalfCol*2.5, height = Width_HalfCol*1.25, pointsize = AxisTxFontSizeSize)
p
dev.off()
openfile(file = graph_path)

# dig sets ####

# set8
wanted1 = c("honeyeaters_pardalote_spatial")
wanted = c("gene", wanted1)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 ) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
# SYN3
set8 = comb_ws2
final_list = tibble('SetNo' = "set8", "gene" = set8)

# set12
wanted1 = c("sunbirds_flowerpecker_gene")
wanted2 = c("sunbirds_flowerpecker_spatial")
wanted = c("gene", wanted1, wanted2)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
# C28H19ORF44, CIB3, EPS15L1, FAM32A, MED26, SLC35E1, TPM4
set12 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set12", "gene" = set12))

# set14
wanted1 = c("hummingbirds_gene")
wanted2 = c("hummingbirds_spatial")
wanted = c("gene", wanted1, wanted2)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
# ERLIN2, EXTL1, GJA9, MYCBP, PAFAH2, PARD3B, RPS6KA2, SMOC2, ZNF703
set14 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set14", "gene" = set14))

# set20
wanted1 = c("honeyeaters_pardalote_gene")
wanted2 = c("honeyeaters_pardalote_spatial")
wanted = c("gene", wanted1, wanted2)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
# ABCB8, ACT85CL, ACTRT2, ATG9B, BPIFCB, C19orf60, COMP, CRLF1, CRTC1, FASTK, KCNH2, KLHL26, KXD1, LOC107055878, LOC107057170, LOC112532681, PRDM16, RTCB, TMEM59L, TMUB1, UPF1
set20 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set20", "gene" = set20))

# set28
wanted1 = c("hummingbirds_absrel")
wanted2 = c("parrots_absrel")
wanted3 = c("honeyeaters_pardalote_absrel")
wanted4 = c("sunbirds_flowerpecker_absrel")
wanted = c("gene", wanted1, wanted2, wanted3, wanted4)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1 &
                                 !!rlang::sym(wanted3) == 1 & !!rlang::sym(wanted4) == 1) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
# CABP1, COL13A1, CTNNAL1, CTSA, CYB561D2, FLVCR2, GRID2IP, HEPACAM2, ICE2, IL1RAP, IRF2, JAM3, KCNK16, KCTD2, MAD1L1, MCOLN2, MLXIPL, MYLK2, MYRF, OAF, OTULIN, PARD6G, PEX11A, PEX6, PRR11, PRUNE1, RAB11FIP4, RALGAPA1, SLC35A2, SLC38A10, SRPK1, STK39, SWT1, TCTE1, TIMM22, TLK2, TMEM8B, TRIM23, TRMT12, UNC13D, VTG1, WBP1L, WRN
set28 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set28", "gene" = set28))


# set15
wanted1 = c("hummingbirds_absrel")
wanted2 = c("parrots_absrel")
wanted3 = c("sunbirds_flowerpecker_absrel")
wanted = c("gene", wanted1, wanted2, wanted3)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1 &
                                 !!rlang::sym(wanted3) == 1 ) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
# AP4S1, ARPC1A, CHST14, GABRA4, IRX6, MALSU1, MDH1, NRF1, NTNG2, PDGFC, RPGR, TNFRSF11A, UBR5
set15 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set15", "gene" = set15))

# set19
wanted1 = c("hummingbirds_absrel")
wanted2 = c("parrots_absrel")
wanted3 = c("honeyeaters_pardalote_absrel")
wanted = c("gene", wanted1, wanted2, wanted3)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1 &
                                 !!rlang::sym(wanted3) == 1 ) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
# AGMO, CD247, COL4A2, DNM3, FAM129B, ILF2, KCNK5, KIF12, LDAH, NHLRC3, NPSR1, OSGIN2, SGCE, SLC4A10, TMEM200B, TNNT3, TTC3, VEGFD, ZNF385B
set19 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set19", "gene" = set19))


# set23
wanted1 = c("hummingbirds_absrel")
wanted2 = c("honeyeaters_pardalote_absrel")
wanted3 = c("sunbirds_flowerpecker_absrel")
wanted = c("gene", wanted1, wanted2, wanted3)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1 &
                                 !!rlang::sym(wanted3) == 1 ) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
# ADAMTS9, CREB3L3, DIS3L2, DNAI1, FAM206A, GAL, GPSM2, HIP1, HNRNPAB, INSYN2, IQCA1, ITGB3, MBD4, MCM2, MKNK1, NUDT15, OMA1, PIK3C2G, POT1, PPP1R21, PPP2R2A, SORCS1, SPATA7, TMX4, TRMT61B, VASN, VPS50
set23 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set23", "gene" = set23))


# set24
wanted1 = c("parrots_absrel")
wanted2 = c("honeyeaters_pardalote_absrel")
wanted3 = c("sunbirds_flowerpecker_absrel")
wanted = c("gene", wanted1, wanted2, wanted3)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1 &
                                 !!rlang::sym(wanted3) == 1 ) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
# ADGRG6, ASB18, C10orf90, CDK15, CNGA2, CREG2, DENND2A, DIABLO, EHD4, EXOSC3, GFRA3, IQSEC1, KIAA1211L, MYBL1, NCKIPSD, PPP1R42, PTPN2, R3HDM1, RAD9B, SCUBE3, SLC26A3, SNX21, TBC1D19, TMEM53, TMPRSS3, UBE3D, WARS2
set24 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set24", "gene" = set24))

# set16
wanted1 = c("parrots_absrel")
wanted2 = c("sunbirds_flowerpecker_absrel")
wanted = c("gene", wanted1, wanted2)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
# C19orf44, CLIC3, DUSP28, EFNB1, GTPBP8, HMGXB3, KIZ, LSM6, MFSD2B, NOV, PGGHG, POLL, RGS16, SFMBT2, SPIK5, TMEM206
set16 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set16", "gene" = set16))

# set21
wanted1 = c("hummingbirds_absrel")
wanted2 = c("parrots_absrel")
wanted = c("gene", wanted1, wanted2)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
# ALDH18A1, DND1, DYTN, FLNB, G6PC, GYG1, HABP2, HAVCR1, HSD3B1, IPO8, IQCJ, MYEF2, PAFAH1B2, PGD, PNPLA3, RGS11, SCAPER, SELENOT, SHOC2, SYDE2, THAP4, TNKS, TRAK1, UBXN11
set21 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set21", "gene" = set21))

# set22
wanted1 = c("parrots_absrel")
wanted2 = c("honeyeaters_pardalote_absrel")
wanted = c("gene", wanted1, wanted2)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
# AASDHPPT, ANK2, ANKK1, ASPA, AVPR1B, BEGAIN, CCDC122, CCDC82, FAM120B, FAM19A4, FBRSL1, FEM1A, FIBCD1, FYB2, GRSF1, IL18RAP, MOB3B, NEBL, PDGFRB, PLEKHD1, RASL11A, RHOT2, RNF220, SDR42E2, TNIP3, TOM1L2
set22 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set22", "gene" = set22))



# set25
wanted1 = c("hummingbirds_absrel")
wanted2 = c("honeyeaters_pardalote_absrel")
wanted = c("gene", wanted1, wanted2)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
set25 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set25", "gene" = set25))


# set26
wanted1 = c("hummingbirds_absrel")
wanted2 = c("sunbirds_flowerpecker_absrel")
wanted = c("gene", wanted1, wanted2)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
set26 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set26", "gene" = set26))

# set27
wanted1 = c("honeyeaters_pardalote_absrel")
wanted2 = c("sunbirds_flowerpecker_absrel")
wanted = c("gene", wanted1, wanted2)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
set27 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set27", "gene" = set27))

# set7
wanted1 = c("hummingbirds_gene")
wanted2 = c("sunbirds_flowerpecker_gene")
wanted = c("gene", wanted1, wanted2)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
# LOC419602
set7 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set7", "gene" = set7))


# set13
wanted1 = c("honeyeaters_pardalote_gene")
wanted2 = c("sunbirds_flowerpecker_gene")
wanted = c("gene", wanted1, wanted2)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
set13 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set13", "gene" = set13))

# set1
wanted1 = c("hummingbirds_absrel")
wanted2 = c("hummingbirds_gene")
wanted3 = c("honeyeaters_pardalote_absrel")
wanted4 = c("sunbirds_flowerpecker_absrel")
wanted = c("gene", wanted1, wanted2, wanted3, wanted4)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1 &
                                 !!rlang::sym(wanted3) == 1 & !!rlang::sym(wanted4) == 1) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
# FANCM
set1 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set1", "gene" = set1))

# set2
wanted1 = c("hummingbirds_spatial")
wanted2 = c("hummingbirds_gene")
wanted3 = c("honeyeaters_pardalote_absrel")
wanted4 = c("sunbirds_flowerpecker_absrel")
wanted = c("gene", wanted1, wanted2, wanted3, wanted4)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1 &
                                 !!rlang::sym(wanted3) == 1 & !!rlang::sym(wanted4) == 1) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
# JADE1
set2 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set2", "gene" = set2))


# set3
wanted1 = c("hummingbirds_spatial")
wanted2 = c("hummingbirds_gene")
wanted3 = c("hummingbirds_absrel")
wanted = c("gene", wanted1, wanted2, wanted3)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1 &
                                 !!rlang::sym(wanted3) == 1 ) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
# SLC30A2
set3 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set3", "gene" = set3))

# set4
wanted1 = c("hummingbirds_absrel")
wanted2 = c("hummingbirds_gene")
wanted3 = c("sunbirds_flowerpecker_absrel")
wanted = c("gene", wanted1, wanted2, wanted3)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1 &
                                 !!rlang::sym(wanted3) == 1 ) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
# RAB11FIP1
set4 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set4", "gene" = set4))

# set5
wanted1 = c("honeyeaters_pardalote_gene")
wanted2 = c("honeyeaters_pardalote_spatial")
wanted3 = c("sunbirds_flowerpecker_gene")
wanted = c("gene", wanted1, wanted2, wanted3)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1 &
                                 !!rlang::sym(wanted3) == 1 ) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
# AOC1
set5 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set5", "gene" = set5))

# set6
wanted1 = c("hummingbirds_absrel")
wanted2 = c("hummingbirds_gene")
wanted = c("gene", wanted1, wanted2)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1 ) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
# AOC1
set6 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set6", "gene" = set6))


# set9
wanted1 = c("hummingbirds_gene")
wanted2 = c("honeyeaters_pardalote_gene")
wanted3 = c("honeyeaters_pardalote_spatial")
wanted = c("gene", wanted1, wanted2, wanted3)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1 &
                                 !!rlang::sym(wanted3) == 1 ) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
# FOXP4L, MDFI
set9 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set9", "gene" = set9))

# set10
wanted1 = c("hummingbirds_gene")
wanted2 = c("honeyeaters_pardalote_absrel")
wanted = c("gene", wanted1, wanted2)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
# NOLC1, TCTE3
set10 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set10", "gene" = set10))

# set11
wanted1 = c("hummingbirds_gene")
wanted2 = c("hummingbirds_spatial")
wanted3 = c("sunbirds_flowerpecker_gene")
wanted = c("gene", wanted1, wanted2, wanted3)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1 & !!rlang::sym(wanted2) == 1 &
                                 !!rlang::sym(wanted3) == 1 ) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
# NOLC1, TCTE3
set11 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set11", "gene" = set11))

# set17
wanted1 = c("honeyeaters_pardalote_gene")
wanted = c("gene", wanted1)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1  ) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
set17 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set17", "gene" = set17))

# set18
wanted1 = c("sunbirds_flowerpecker_gene")
wanted = c("gene", wanted1)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1  ) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
set18 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set18", "gene" = set18))

# set29
wanted1 = c("hummingbirds_gene")
wanted = c("gene", wanted1)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1  ) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
set29 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set29", "gene" = set29))

# set30
wanted1 = c("parrots_absrel")
wanted = c("gene", wanted1)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1  ) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
set30 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set30", "gene" = set30))

# set31
wanted1 = c("honeyeaters_pardalote_absrel")
wanted = c("gene", wanted1)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1  ) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
set31 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set31", "gene" = set31))

# set32
wanted1 = c("sunbirds_flowerpecker_absrel")
wanted = c("gene", wanted1)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1  ) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
set32 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set32", "gene" = set32))

# set33
wanted1 = c("hummingbirds_absrel")
wanted = c("gene", wanted1)
other = setdiff(names(comb_wide), wanted)
comb_ws = comb_wide %>% filter(!!rlang::sym(wanted1) == 1  ) %>%
  dplyr::select(all_of(wanted), all_of(other))
sumt = apply(comb_ws[(length(wanted)+1):ncol(comb_ws)], 1, sum)
comb_ws2 = comb_ws %>% 
  mutate(sumt = sumt) %>%
  filter(sumt == 0) %>%
  arrange(gene) %>% pull(gene)
length(comb_ws2)
paste0(comb_ws2, collapse = ", ")
cat(paste0(comb_ws2, collapse = "\n"))
set33 = comb_ws2
final_list = final_list %>% bind_rows(tibble('SetNo' = "set33", "gene" = set33))


# check
length(unique(final_list$SetNo))
setdiff(paste0("set", 1:33), unique(final_list$SetNo))
nrow(final_list) == length(unique(comb$gene))

final_list2 = final_list %>%
  mutate(set_notes = ifelse(SetNo == "set6", "overlap - 2 analyses (per clade)", NA),
         set_notes = ifelse(SetNo == "set3", "overlap - 3 analyses (per clade)", set_notes),
         set_notes = ifelse(SetNo == "set8", "no overlap - honeyeaters_pardalote_spatial", 
                            set_notes),
         set_notes = ifelse(SetNo %in% paste0("set", c(12, 14, 20)), 
                            "overlap - 2 analyses (per clade)", 
                            set_notes),
         set_notes = ifelse(SetNo %in% paste0("set", c(28)), 
                            "absrel overlap - 4 clades", 
                            set_notes),
         set_notes = ifelse(SetNo %in% paste0("set", c(15, 19, 23, 24)), 
                            "absrel overlap - 3 clades", 
                            set_notes),
         set_notes = ifelse(SetNo %in% paste0("set", c(16, 21, 22, 25, 26, 27)), 
                            "absrel overlap - 2 clades", 
                            set_notes),
         set_notes = ifelse(SetNo %in% paste0("set", c(7, 13)), 
                            "gene overlap - 2 clades", 
                            set_notes),
         
         set_notes = ifelse(SetNo == "set29", "no overlap - hummingbirds_gene", 
                           set_notes),
         set_notes = ifelse(SetNo == "set30", "no overlap - parrots_absrel", 
                            set_notes),
         set_notes = ifelse(SetNo == "set31", "no overlap - honeyeaters_pardalote_absrel", 
                            set_notes),
         set_notes = ifelse(SetNo == "set32", "no overlap - sunbirds_flowerpecker_absrel", 
                            set_notes),
         set_notes = ifelse(SetNo == "set33", "no overlap - hummingbirds_absrel", 
                            set_notes),
         set_notes = ifelse(SetNo == "set17", "no overlap - honeyeaters_pardalote_gene", 
                           set_notes),
         set_notes = ifelse(SetNo == "set18", "no overlap - sunbirds_flowerpecker_gene", 
                            set_notes),
         
         set_notes = ifelse(SetNo %in% paste0("set", c(9, 11)), 
                            "gene overlap - 2 clades & spatial", 
                            set_notes),
         set_notes = ifelse(SetNo %in% paste0("set", c(1)), 
                            "overlap - 2 analyses, 3 clades", 
                            set_notes),
         set_notes = ifelse(SetNo %in% paste0("set", c(2)), 
                            "overlap - 3 analyses, 3 clades", 
                            set_notes),
         set_notes = ifelse(SetNo %in% paste0("set", c(4, 5)), 
                            "overlap - 2 analyses, 2 clades", 
                            set_notes),
         set_notes = ifelse(SetNo %in% paste0("set", c(10)), 
                            "overlap - 2 analyses, 2 clades", 
                            set_notes)
         )
final_list3 = final_list2 %>% filter(is.na(set_notes)) %>% distinct(SetNo); final_list3


# get anno #### 
ensembl = useMart("ensembl")
# datasets = listDatasets(ensembl)
# head(datasets) 
# dim(datasets) # 203 3
# datasets$dataset[grep("gal", datasets$dataset, ignore.case = T)]

ensembl_gg = useDataset("ggallus_gene_ensembl", mart=ensembl)
filters = listFilters(ensembl_gg) # query # from
filters$name[grep('symbo', filters$name)] # hgnc_symbol

attributes = listAttributes(ensembl_gg) # to
attributes$name[grep('name', attributes$name)] # wikigene_name external_gene_name
attributes$name[grep('ggallus', attributes$name)] # 

geneDescr =  getBM(attributes=c('hgnc_symbol', 'description'), 
                   filters = 'hgnc_symbol', 
                   values = final_list2$gene, 
                   mart = ensembl_gg)
geneDescr2 = geneDescr %>%
  separate(description, into = c("description", "source"), sep = " \\[Source:", remove = F) %>%
  separate(source, into = c("source", "acc"), sep = ";Acc:", remove = T) %>%
  mutate(acc = gsub("\\]", "", acc),
         acc = gsub("HGNC:", "", acc))

go_geneDescr = final_list2 %>% 
  dplyr::left_join(geneDescr2, by = c("gene" = "hgnc_symbol")) %>%
  arrange(SetNo, gene) %>%
  dplyr::select(SetNo, set_notes, gene, description, source, acc)
nrow(go_geneDescr) # 1013
good = go_geneDescr %>% filter(! is.na(description))


# look up in human
tolookagain = go_geneDescr %>% filter(is.na(description))
datasets = listDatasets(ensembl)
head(datasets)
dim(datasets) # 203 3
datasets$dataset[grep("hsa", datasets$dataset, ignore.case = T)]

ensembl_h = useDataset("hsapiens_gene_ensembl", mart=ensembl)
filters = listFilters(ensembl_h) # query # from
filters$name[grep('symbo', filters$name)] # hgnc_symbol

attributes = listAttributes(ensembl_h) # to
attributes$name[grep('name', attributes$name)] # wikigene_name external_gene_name
attributes$name[grep('hsa', attributes$name)] # 

geneDescr =  getBM(attributes=c('hgnc_symbol', 'description'), 
                   filters = 'hgnc_symbol', 
                   values = tolookagain$gene, 
                   mart = ensembl_h)
geneDescr2 = geneDescr %>%
  separate(description, into = c("description", "source"), sep = " \\[Source:", remove = F) %>%
  separate(source, into = c("source", "acc"), sep = ";Acc:", remove = T) %>%
  mutate(acc = gsub("\\]", "", acc),
         acc = gsub("HGNC:", "", acc))

go_geneDescr2 = tolookagain %>% dplyr::select(SetNo, set_notes, gene) %>% 
  left_join(geneDescr2, by = c("gene" = "hgnc_symbol")) 

good2 = good %>% bind_rows(go_geneDescr2) %>%
  mutate(SetNo = gsub("set", "", SetNo), 
         SetNo = as.numeric(SetNo)) %>% 
  arrange(SetNo, gene)
nrow(good2)

path = paste0(project_path, "lists_absrel_gene_spatial.tsv")
write_tsv(good2, path, na = "")
