
# setup ouch
# .libPaths()
# libpath = .libPaths()[2]
# install.packages("devtools", lib = libpath)
# library(devtools)
# install_github("kingaa/ouch", dependencies = T)
# install_github("kingaa/ouch", lib = libpath)
# BiocManager::install("edgeR", lib = libpath)
# library(edgeR) # for normalization


library(tidyverse)
library(ouch)
library(parallel)

# args <- commandArgs(trailingOnly = TRUE)
# #args 1 is the tissue
# #args 2 is number of cores
# tmp_tissue = args[1]
# core = args[2] 
# tmp_tissue = "duodenum"
core = 12

input_path = paste0(getwd(), "/Seq_Data/globus/rnaseq/00_inputs/")
project_path = paste0(getwd(), "/Seq_Data/globus/rnaseq/01_ouch/")

# input_path = paste0("/home/mpg08/mko/Nectar/analysis/rnaseq/00_inputs/")
# project_path = paste0("/home/mpg08/mko/Nectar/analysis/rnaseq/01_ouch/")

# input = c("duodenum", "heart","pectoralis" ,"liver") 
input = c("palate") 

sapply(input, function(tis){

  tmp_tissue = tis
  output_path = paste0(project_path, "/output/", tmp_tissue, "/")
  dir.create(output_path, recursive = T)
  
  # 1. load data ====
  # path = paste0(project_path, "NectarTree02.nw") # no branch lengths # 
  path = paste0(input_path, "nonconserved_4d.nw") 
  tree = ape::read.tree(path)
  
  tree$tip.label = gsub("HLcalAnn5", "Anna's_hummingbird", tree$tip.label)
  tree$tip.label = gsub("HLnymHol2", "cockatiel", tree$tip.label)
  tree$tip.label = gsub("HLapuApu1", "common_swift", tree$tip.label)
  tree$tip.label = gsub("HLphyNov1", "New_Holland_honeyeater", tree$tip.label)
  tree$tip.label = gsub("HLtriMol2", "rainbow_lorikeet", tree$tip.label)
  tree$tip.label = gsub("HLtaeGut4", "zebra_finch", tree$tip.label)

  
  # 2. ouch ====
  set.seed(100)
  
  path = list.files(input_path, paste0("average_gene_expr_TMM_", tmp_tissue), full.names = T)
  tmp_exp_sel = read_tsv(path)
  print(tmp_exp_sel$bird)
  
  tmp_matrix = tmp_exp_sel[, -1] %>% as.matrix
  row.names(tmp_matrix) = tmp_exp_sel$bird
  data_mean.norm_transp = tmp_matrix
  
  ### make an ouchtree out of the phy-format tree
  tree_m = tree %>% ape::keep.tip(tmp_exp_sel$bird) # keep tips with rnaseq data
  ot <- ape2ouch(tree_m)

  run_ouch = function(iCol){
    # iCol = 409
    # print(iCol)
    tmp_gene = colnames(data_mean.norm_transp)[iCol]
    tmp_gene_matrix = data_mean.norm_transp[, iCol] %>% as.matrix(ncol = 1)
    colnames(tmp_gene_matrix) = tmp_gene
    
    if ( length(unique(unique(tmp_gene_matrix))) >1 ){
      
      ### merge data with tree info
      otd <- as(ot,"data.frame")
      otd <- merge(otd, tmp_gene_matrix, by.x = "labels",by.y = "row.names",all = TRUE)
      ### row-names are used by 'hansen'
      rownames(otd) <- otd$nodes
      # print(otd)
      ### this data-frame now contains the data as well as the tree geometry
      
      ### now remake the ouch tree
      ot2 <- with(otd, ouchtree(nodes = nodes,ancestors = ancestors,times = times,labels = labels))
      # plot(ot2, node.names = TRUE)
      
      # brownian
      b1 <- brown(tree = ot2, data = otd[tmp_gene])
      # summary(b1)
      
      ### evaluate an OU model with a single, global selective regime
      otd$regimes <- as.factor("global")
      hg <- try(hansen(
        tree = ot2,
        data = otd[tmp_gene],
        regimes = otd["regimes"],
        sqrt.alpha =  1,
        sigma = 1,
        # maxit=10000
        reltol=1e-4
      ), TRUE)
      # summary(hg)
      # plot(hg)
      
      ### evaluate an OU model with a selective regime on nectar species
      ### paint regime
      r <- paint(ot2, branch = c("7" = "nectar", "9" = "nectar", "10" = "nectar"))
      # plot(ot2, regimes = r, node.names = TRUE)
      
      otd$regimes <- r
      hn <- try(hansen(
        tree = ot2,
        data = otd[tmp_gene],
        regimes = otd["regimes"],
        sqrt.alpha = 1,
        sigma = 1,
        reltol=1e-4
      ), TRUE)
      # summary(hn)
      # plot(hn)
      
      if(  !inherits(hg, "try-error") & !inherits(hn, "try-error") ){
        sigma_Br = coef(b1)$sigma.sq.matrix[,1] # sigma.squared
        theta_Br = coef(b1)$theta[[1]]
        logLik_Br = logLik(b1)
        aic_Br = summary(b1)$aic
        aicc_Br = summary(b1)$aic.c # small sample bias-correction
        bic_Br = summary(b1)$sic
        
        sigma_Hg = coef(hg)$sigma.sq.matrix[,1]
        alpha_Hg = coef(hg)$alpha.matrix[,1]
        theta_Hg = coef(hg)$theta[[1]]
        logLik_Hg = logLik(hg)
        pvalue_Hg =  1 - pchisq((logLik(hg) - logLik(b1)) * 2, 1)
        aic_Hg = summary(hg)$aic
        aicc_Hg = summary(hg)$aic.c # small sample bias-correction
        bic_Hg = summary(hg)$sic
        
        sigma_Hn = coef(hn)$sigma.sq.matrix[,1]
        alpha_Hn = coef(hn)$alpha.matrix[,1]
        theta_Hn_nectar = coef(hn)$theta[[1]][1]
        theta_Hn_nonnectar = coef(hn)$theta[[1]][2]
        logLik_Hn = logLik(hn)
        pvalue_Hn =  1 - pchisq((logLik(hn) - logLik(b1)) * 2, summary(hn)$dof-1)
        aic_Hn = summary(hn)$aic
        aicc_Hn = summary(hn)$aic.c # small sample bias-correction
        bic_Hn = summary(hn)$sic
        
        pvalue_HnHg =  1 - pchisq((logLik(hn) - logLik(hg)) * 2, summary(hn)$dof-1)
        # k <- summary(hn)$dof
        # n = nrow(tmp_gene_matrix)
        # aicc = -2*logLik_Hn + 2*k*n/(n-k-1) # this equals to aicc_Hn
        
        tb = tibble("gene" = tmp_gene, 
                    "sigma_Br" = sigma_Br, "theta_Br" = theta_Br, "logLik_Br" = logLik_Br, 
                    "aic_Br" = aic_Br, "aicc_Br" = aicc_Br, "bic_Br" = bic_Br,
                    "sigma_Hg" = sigma_Hg, "alpha_Hg" = alpha_Hg, "theta_Hg" = theta_Hg, 
                    "logLik_Hg" = logLik_Hg, "pvalue_Hg" = pvalue_Hg,
                    "aic_Hg" = aic_Hg, "aicc_Hg" = aicc_Hg, "bic_Hg" = bic_Hg,
                    "sigma_Hn" = sigma_Hn, "alpha_Hn" = alpha_Hn, 
                    "theta_Hn_nectar" = theta_Hn_nectar, "theta_Hn_nonnectar" = theta_Hn_nonnectar,
                    "logLik_Hn" = logLik_Hn, "pvalue_Hn" = pvalue_Hn, 
                    "aic_Hn" = aic_Hn, "aicc_Hn" = aicc_Hn, "bic_Hn" = bic_Hn,
                    "pvalue_HnHg" = pvalue_HnHg) 
        return(tb)
      }
    }#else{print(iCol)}
  }
  
  i_to_try = 1:ncol(tmp_matrix) # iCol of data_mean.norm_transp
  
  # Start the clock!
  ptm1 <- proc.time()
  
  # out_ouch = mapply(run_ouch, iCol = i_to_try, SIMPLIFY = F)
  out_ouch = mclapply(X = i_to_try, FUN = run_ouch, mc.cores= core, mc.preschedule = FALSE)
  
  # Stop the clock
  pt1 = proc.time() - ptm1; pt1
  
  out_ouch_tb = out_ouch %>% bind_rows() %>% 
    mutate(qvalue_Hg = p.adjust(pvalue_Hg, method="fdr"),
           qvalue_Hn = p.adjust(pvalue_Hn, method="fdr"),
           qvalue_HnHg = p.adjust(pvalue_HnHg, method="fdr"),
           var_Hg = sigma_Hg / (alpha_Hg*2),
           var_Hn = sigma_Hn / (alpha_Hn*2),
           tissue = tmp_tissue
    ) %>% 
    dplyr::select(tissue, gene, matches("Br"), matches("Hg"), matches("Hn"))
  # dim(out_ouch_tb)
  out_ouch_tbm = out_ouch_tb %>% as.matrix(nrow = length(i_to_try)) %>% as_tibble()
  # dim(out_ouch_tbm)
  
  path = paste0(output_path, tmp_tissue, "_ouch_TMM_", Sys.Date(), ".tsv")
  write_tsv( out_ouch_tbm, path, na = "")

})

