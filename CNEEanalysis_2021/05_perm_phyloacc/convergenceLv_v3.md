4 binomial distributions
================
maggie
2022-11-05

-   <a href="#all-cnees-within-100kb-of-a-gene"
    id="toc-all-cnees-within-100kb-of-a-gene">all cnees within 100kb of a
    gene</a>
-   <a href="#kegg-annotation" id="toc-kegg-annotation">kegg annotation</a>
-   <a href="#nectar-pathway-gene-and-acc-cnees"
    id="toc-nectar-pathway-gene-and-acc-cnees">nectar: pathway, gene, and
    acc cnees</a>
-   <a href="#core-non-nectar-pathway-gene-and-acc-cnees"
    id="toc-core-non-nectar-pathway-gene-and-acc-cnees">core non-nectar:
    pathway, gene, and acc cnees</a>
-   <a href="#probability" id="toc-probability">Probability</a>
-   <a href="#nectar-binomial" id="toc-nectar-binomial">nectar: Binomial</a>
-   <a href="#core-non-nectar-binomial"
    id="toc-core-non-nectar-binomial">core non-nectar: Binomial</a>
-   <a href="#comparison" id="toc-comparison">comparison</a>

## all cnees within 100kb of a gene

``` r
set.seed(100)
library(tidyverse) 
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ## ✔ ggplot2 3.3.6      ✔ purrr   0.3.5 
    ## ✔ tibble  3.1.8      ✔ dplyr   1.0.10
    ## ✔ tidyr   1.2.1      ✔ stringr 1.4.1 
    ## ✔ readr   2.1.3      ✔ forcats 0.5.2 
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(ggvenn)
```

    ## Loading required package: grid

``` r
project_path = paste0(getwd(), "/postPhyloAcc/")
out_path = paste0(project_path, "convergentLv/")

## all cnees within 100kb of a gene (created using bedtools)
path = paste0(project_path, "linear/cnee_ncbigene100kb_intersect_allcnees.bed")
allcnees = read_tsv(path, col_names = F) %>% 
  filter(X6 == "protein_coding") %>%
  dplyr::select(X4, X5, X10) %>%
  dplyr::rename("gene" = "X4", "id" = "X10", "combo" = "X5") %>%
  separate(combo, into = c(NA, "ncbi"), sep = "GeneID:", remove = F) %>% 
  dplyr::select(gene, id) 
```

    ## Rows: 992273 Columns: 10
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (6): X1, X4, X5, X6, X7, X10
    ## dbl (4): X2, X3, X8, X9
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
dim(allcnees) # 613780      2
```

    ## [1] 613780      2

``` r
head(allcnees)
```

    ## # A tibble: 6 × 2
    ##   gene      id   
    ##   <chr>     <chr>
    ## 1 LOC430443 CNEE3
    ## 2 LOC430443 CNEE4
    ## 3 LOC430443 CNEE5
    ## 4 GOLGB1    CNEE3
    ## 5 GOLGB1    CNEE4
    ## 6 GOLGB1    CNEE5

``` r
# gene      id   
# <chr>     <chr>
# 1 LOC430443 CNEE3
```

## kegg annotation

``` r
## custom kegg annotation
path = paste0(project_path, "ncbi_anno_kegg_hsa_2022-10-22.tsv")
custom_anno = read_tsv(path)
```

    ## Rows: 17340 Columns: 7
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (5): gene, ensembl_gene_id, hsapiens_homolog_ensembl_gene, pathway_hsa, ...
    ## dbl (2): entrezgene_id_gal, entrez_hsa
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
custom_anno_allcnees = custom_anno %>% left_join(allcnees) %>% 
  distinct()
```

    ## Joining, by = "gene"

``` r
dim(custom_anno_allcnees) # 819389      8
```

    ## [1] 819389      8

``` r
head(custom_anno_allcnees)
```

    ## # A tibble: 6 × 8
    ##   gene      entrezgene_id_gal ensembl_ge…¹ hsapi…² entre…³ pathw…⁴ descr…⁵ id   
    ##   <chr>                 <dbl> <chr>        <chr>     <dbl> <chr>   <chr>   <chr>
    ## 1 LAMTOR1           100857186 ENSGALG0000… ENSG00…   55004 path:h… mTOR s… CNEE…
    ## 2 LAMTOR1           100857186 ENSGALG0000… ENSG00…   55004 path:h… mTOR s… CNEE…
    ## 3 LAMTOR1           100857186 ENSGALG0000… ENSG00…   55004 path:h… mTOR s… CNEE…
    ## 4 RAB11FIP3         100857269 ENSGALG0000… ENSG00…    9727 path:h… Endocy… CNEE…
    ## 5 RAB11FIP3         100857269 ENSGALG0000… ENSG00…    9727 path:h… Endocy… CNEE…
    ## 6 RAB11FIP3         100857269 ENSGALG0000… ENSG00…    9727 path:h… Endocy… CNEE…
    ## # … with abbreviated variable names ¹​ensembl_gene_id,
    ## #   ²​hsapiens_homolog_ensembl_gene, ³​entrez_hsa, ⁴​pathway_hsa, ⁵​description_hsa

``` r
## kegg pathway id and description
custom_anno_id = custom_anno %>% 
  dplyr::select(pathway_hsa, description_hsa) %>% distinct()
head(custom_anno_id)
```

    ## # A tibble: 6 × 2
    ##   pathway_hsa   description_hsa                      
    ##   <chr>         <chr>                                
    ## 1 path:hsa04150 mTOR signaling pathway               
    ## 2 path:hsa04144 Endocytosis                          
    ## 3 path:hsa04070 Phosphatidylinositol signaling system
    ## 4 path:hsa04110 Cell cycle                           
    ## 5 path:hsa04114 Oocyte meiosis                       
    ## 6 path:hsa04120 Ubiquitin mediated proteolysis

``` r
## total cnees of a kegg pathway
custom_anno_allcnees_sum = custom_anno_allcnees %>% 
  group_by(pathway_hsa) %>% 
  summarise(total_cnees = length(unique(id)))
dim(custom_anno_allcnees_sum) # 350 2
```

    ## [1] 350   2

``` r
head(custom_anno_allcnees_sum)
```

    ## # A tibble: 6 × 2
    ##   pathway_hsa   total_cnees
    ##   <chr>               <int>
    ## 1 path:hsa00010         689
    ## 2 path:hsa00020         121
    ## 3 path:hsa00030         459
    ## 4 path:hsa00040         846
    ## 5 path:hsa00051         957
    ## 6 path:hsa00052         777

## nectar: pathway, gene, and acc cnees

``` r
## nectar: pathway, gene, and acc cnees 
path = paste0(out_path, "convergent_pathwayLv_full_2022-10-24.tsv")
pathway_acc_anno = read_tsv(path)
```

    ## Rows: 7338 Columns: 12
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (7): group, gene, id, ensembl_gene_id, hsapiens_homolog_ensembl_gene, pa...
    ## dbl (5): ncbi, accel, total, entrezgene_id_gal, entrez_hsa
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
## number of acc cnees for each pathway for each clade
pathway_acc_sum = pathway_acc_anno %>% group_by(pathway_hsa, description_hsa, group) %>% 
  summarise(total_acccnees = length(unique(id))) %>% 
  left_join(custom_anno_allcnees_sum)
```

    ## `summarise()` has grouped output by 'pathway_hsa', 'description_hsa'. You can
    ## override using the `.groups` argument.
    ## Joining, by = "pathway_hsa"

``` r
head(pathway_acc_sum)
```

    ## # A tibble: 6 × 5
    ## # Groups:   pathway_hsa, description_hsa [2]
    ##   pathway_hsa   description_hsa              group               total…¹ total…²
    ##   <chr>         <chr>                        <chr>                 <int>   <int>
    ## 1 path:hsa00010 Glycolysis / Gluconeogenesis honeyeaters               2     689
    ## 2 path:hsa00010 Glycolysis / Gluconeogenesis hummingbirds              7     689
    ## 3 path:hsa00010 Glycolysis / Gluconeogenesis parrots                   2     689
    ## 4 path:hsa00010 Glycolysis / Gluconeogenesis sunbirds_flowerpec…       1     689
    ## 5 path:hsa00030 Pentose phosphate pathway    hummingbirds              2     459
    ## 6 path:hsa00030 Pentose phosphate pathway    parrots                   1     459
    ## # … with abbreviated variable names ¹​total_acccnees, ²​total_cnees

``` r
## total cnees per pathway
pathway_acc_sum2 = pathway_acc_sum %>% 
  dplyr::select(pathway_hsa, total_cnees) %>% distinct()
```

    ## Adding missing grouping variables: `description_hsa`

``` r
head(pathway_acc_sum2)
```

    ## # A tibble: 6 × 3
    ## # Groups:   pathway_hsa, description_hsa [6]
    ##   description_hsa                          pathway_hsa   total_cnees
    ##   <chr>                                    <chr>               <int>
    ## 1 Glycolysis / Gluconeogenesis             path:hsa00010         689
    ## 2 Pentose phosphate pathway                path:hsa00030         459
    ## 3 Pentose and glucuronate interconversions path:hsa00040         846
    ## 4 Fructose and mannose metabolism          path:hsa00051         957
    ## 5 Galactose metabolism                     path:hsa00052         777
    ## 6 Ascorbate and aldarate metabolism        path:hsa00053         358

``` r
## summarizing convergent pathway
pathwaylv_conv = pathway_acc_anno %>% 
  # keep only pathway lv
  # dplyr::select(group, pathway_hsa, description_hsa) %>% distinct() %>% 
  group_by(pathway_hsa) %>% 
  mutate(n_clade = length(unique(group))) %>% 
  dplyr::select(pathway_hsa, description_hsa, n_clade) %>% 
  distinct()
table(pathwaylv_conv$n_clade)
```

    ## 
    ##   1   2   3   4 
    ##  14  42  55 225

``` r
# 1   2   3   4 
# 14  42  55 225 
head(pathwaylv_conv)
```

    ## # A tibble: 6 × 3
    ## # Groups:   pathway_hsa [6]
    ##   pathway_hsa   description_hsa               n_clade
    ##   <chr>         <chr>                           <int>
    ## 1 path:hsa00500 Starch and sucrose metabolism       4
    ## 2 path:hsa01100 Metabolic pathways                  4
    ## 3 path:hsa04151 PI3K-Akt signaling pathway          4
    ## 4 path:hsa04152 AMPK signaling pathway              4
    ## 5 path:hsa04910 Insulin signaling pathway           4
    ## 6 path:hsa04922 Glucagon signaling pathway          4

## core non-nectar: pathway, gene, and acc cnees

``` r
## core non-nectar: pathway, gene, and acc cnees 
path = paste0(out_path, "convergent_pathwayLv_control_2022-11-05.tsv")
pathway_acc_anno_con = read_tsv(path)
```

    ## Rows: 17148 Columns: 12
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (7): group, gene, id, ensembl_gene_id, hsapiens_homolog_ensembl_gene, pa...
    ## dbl (5): ncbi, accel, total, entrezgene_id_gal, entrez_hsa
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
## number of acc cnees for each pathway for each clade
pathway_acc_con_sum = pathway_acc_anno_con %>% group_by(pathway_hsa, description_hsa, group) %>% 
  summarise(total_acccnees = length(unique(id))) %>% 
  left_join(custom_anno_allcnees_sum)
```

    ## `summarise()` has grouped output by 'pathway_hsa', 'description_hsa'. You can
    ## override using the `.groups` argument.
    ## Joining, by = "pathway_hsa"

``` r
head(pathway_acc_con_sum)
```

    ## # A tibble: 6 × 5
    ## # Groups:   pathway_hsa, description_hsa [2]
    ##   pathway_hsa   description_hsa              group      total_acccnees total_c…¹
    ##   <chr>         <chr>                        <chr>               <int>     <int>
    ## 1 path:hsa00010 Glycolysis / Gluconeogenesis falcons                 3       689
    ## 2 path:hsa00010 Glycolysis / Gluconeogenesis lyrebirds               4       689
    ## 3 path:hsa00010 Glycolysis / Gluconeogenesis passerides              3       689
    ## 4 path:hsa00010 Glycolysis / Gluconeogenesis swifts                  8       689
    ## 5 path:hsa00020 Citrate cycle (TCA cycle)    lyrebirds               1       121
    ## 6 path:hsa00020 Citrate cycle (TCA cycle)    passerides              1       121
    ## # … with abbreviated variable name ¹​total_cnees

``` r
## total cnees per pathway
pathway_acc_con_sum2 = pathway_acc_con_sum %>% 
  dplyr::select(pathway_hsa, total_cnees) %>% distinct()
```

    ## Adding missing grouping variables: `description_hsa`

``` r
head(pathway_acc_con_sum2)
```

    ## # A tibble: 6 × 3
    ## # Groups:   pathway_hsa, description_hsa [6]
    ##   description_hsa                          pathway_hsa   total_cnees
    ##   <chr>                                    <chr>               <int>
    ## 1 Glycolysis / Gluconeogenesis             path:hsa00010         689
    ## 2 Citrate cycle (TCA cycle)                path:hsa00020         121
    ## 3 Pentose phosphate pathway                path:hsa00030         459
    ## 4 Pentose and glucuronate interconversions path:hsa00040         846
    ## 5 Fructose and mannose metabolism          path:hsa00051         957
    ## 6 Galactose metabolism                     path:hsa00052         777

``` r
## summarizing convergent pathway
pathwaylv_conv_con = pathway_acc_anno_con %>% 
  # keep only pathway lv
  # dplyr::select(group, pathway_hsa, description_hsa) %>% distinct() %>% 
  group_by(pathway_hsa) %>% 
  mutate(n_clade = length(unique(group))) %>% 
  dplyr::select(pathway_hsa, description_hsa, n_clade) %>% 
  distinct()
table(pathwaylv_conv_con$n_clade)
```

    ## 
    ##   1   2   3   4 
    ##   6  13  36 282

``` r
  # 1   2   3   4 
  # 6  13  36 282
head(pathwaylv_conv_con)
```

    ## # A tibble: 6 × 3
    ## # Groups:   pathway_hsa [6]
    ##   pathway_hsa   description_hsa                         n_clade
    ##   <chr>         <chr>                                     <int>
    ## 1 path:hsa04072 Phospholipase D signaling pathway             4
    ## 2 path:hsa04080 Neuroactive ligand-receptor interaction       4
    ## 3 path:hsa04724 Glutamatergic synapse                         4
    ## 4 path:hsa05030 Cocaine addiction                             4
    ## 5 path:hsa04910 Insulin signaling pathway                     4
    ## 6 path:hsa03010 Ribosome                                      4

## Probability

``` r
# num of total cnees 
t = 363747 
# num of cnees accelerated in any branch (logBF3 >=10)
# t = 61468 # same results

## probability
# nectar
hump = 1691/t
parp = 1125/t
honp = 523/t
sunp = 481/t

# non-nectar
p = 2418/t # passerides 
l = 1591/t # lyrebirds
f = 1863/t # falcons
sw = 2310/t # swifts
```

## nectar: Binomial

``` r
tmp_pathways = unique(pathway_acc_sum$pathway_hsa)
length(tmp_pathways) # 336
```

    ## [1] 336

``` r
sp_seq = c("hummingbirds", "parrots", "honeyeaters", "sunbirds_flowerpecker") # as rmultinom
N_random_samples = 100000

x_nectar = sapply(tmp_pathways, function(s){
  # s = tmp_pathways[2]

  # for each pathway, get the number of acccnee for each clade
  tmp = pathway_acc_sum %>% filter(pathway_hsa == s)
  tmp_matrix = tmp[match(sp_seq, tmp$group), ] # sort clades (row) as rmultinom
  tmp_matrix_acc = tmp_matrix$total_acccnees
  tmp_matrix_acc[is.na(tmp_matrix_acc)] <- 0   # for clades without any acc. cnee, replace NA to 0
  
  # random sample
  N_cnee_pathway = unique(tmp$total_cnees)
  # print(N_cnee_pathway)
  
  # number of observations, number of trials (zero or more).
  rb_hum = rbinom(N_random_samples, N_cnee_pathway, hump) 
  rb_par = rbinom(N_random_samples, N_cnee_pathway, parp)
  rb_honey = rbinom(N_random_samples, N_cnee_pathway, honp)
  rb_sun = rbinom(N_random_samples, N_cnee_pathway, sunp)
  nectar_bi = rbind(rb_hum, rb_par, rb_honey, rb_sun)
  # dim(nectar_bi) # 4  100
  # nectar_bi[, 1:7]
  
  # for each random sample (a col in nectar_md), 
  # check if the expected is greater than or equal to the observed (num. acc. cnees) for each clade (row)
  ind <- ( nectar_bi >= tmp_matrix_acc ) 
  sum_col = apply(ind, 2, sum) # sum up by col, 
  # if all 4 clades' expected is >= observed, the col will have a sum of 4
  r = sum(sum_col == 4) # sum up how many random samples have more greater than or equal to the observed
  return(r)

})

x_nectar_tb = x_nectar %>% as_tibble(rownames = "pathway_hsa") %>% 
  mutate(p_val = (value+1)/(N_random_samples+1)) %>%  #  (r+1)/(n+1)
  left_join(custom_anno_id) %>% 
  left_join(pathway_acc_con_sum2)
```

    ## Joining, by = "pathway_hsa"
    ## Joining, by = c("pathway_hsa", "description_hsa")

``` r
head(x_nectar_tb)
```

    ## # A tibble: 6 × 5
    ##   pathway_hsa   value   p_val description_hsa                          total_c…¹
    ##   <chr>         <int>   <dbl> <chr>                                        <int>
    ## 1 path:hsa00010   433 0.00434 Glycolysis / Gluconeogenesis                   689
    ## 2 path:hsa00030  5976 0.0598  Pentose phosphate pathway                      459
    ## 3 path:hsa00040 51877 0.519   Pentose and glucuronate interconversions       846
    ## 4 path:hsa00051   694 0.00695 Fructose and mannose metabolism                957
    ## 5 path:hsa00052   956 0.00957 Galactose metabolism                           777
    ## 6 path:hsa00053  5930 0.0593  Ascorbate and aldarate metabolism              358
    ## # … with abbreviated variable name ¹​total_cnees

``` r
hist(x_nectar_tb$p_val)
```

![](convergenceLv_v3_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
hist(x_nectar_tb$total_cnees)
```

![](convergenceLv_v3_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
ggplot(x_nectar_tb, aes(x = total_cnees, y = p_val)) +
  geom_point() 
```

    ## Warning: Removed 5 rows containing missing values (geom_point).

![](convergenceLv_v3_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

``` r
ggplot(x_nectar_tb, aes(x = total_cnees, y = p_val)) +
  geom_point() +
  scale_x_continuous(limits = c(NA, 2500))
```

    ## Warning: Removed 100 rows containing missing values (geom_point).

![](convergenceLv_v3_files/figure-gfm/unnamed-chunk-6-4.png)<!-- -->

``` r
x_nectar_tb_sig = x_nectar_tb %>% 
  filter(p_val<=0.05)
nrow(x_nectar_tb_sig) # p_val<= 0.05: 153; p_val<= 0.01: 79
```

    ## [1] 153

## core non-nectar: Binomial

``` r
tmp_pathways = unique(pathway_acc_con_sum$pathway_hsa)
length(tmp_pathways) # 337
```

    ## [1] 337

``` r
sp_seq = c("swifts", "falcons", "lyrebirds", "passerides")
N_random_samples = 100000

x_nonnectar = sapply(tmp_pathways, function(s){
  # s = tmp_pathways[2]
  
   # for each pathway, get the number of acccnee for each clade
  tmp = pathway_acc_con_sum %>% filter(pathway_hsa == s)
  tmp_matrix = tmp[match(sp_seq, tmp$group), ] # sort clades (row) as rmultinom
  tmp_matrix_acc = tmp_matrix$total_acccnees
  tmp_matrix_acc[is.na(tmp_matrix_acc)] <- 0   # for clades without any acc. cnee, replace NA to 0
  
  # random sample
  N_cnee_pathway = unique(tmp$total_cnees)
  # print(N_cnee_pathway)
  
  # number of observations, number of trials (zero or more).
  rb_passerides = rbinom(N_random_samples, N_cnee_pathway, p) 
  rb_lyrebirds = rbinom(N_random_samples, N_cnee_pathway, l)
  rb_falcons = rbinom(N_random_samples, N_cnee_pathway, f)
  rb_swifts = rbinom(N_random_samples, N_cnee_pathway, sw)
  non_nectar_bi = rbind(rb_passerides, rb_lyrebirds, rb_falcons, rb_swifts) 

  # for each random sample (a col in nectar_md), 
  # check if the expected is greater than or equal to the observed (num. acc. cnees) for each clade (row)
  ind <- ( non_nectar_bi >= tmp_matrix_acc ) 
  sum_col = apply(ind, 2, sum) # sum up by col, 
  # if all 4 clades' expected is >= observed, the col will have a sum of 4
  r = sum(sum_col == 4) # sum up how many random samples have more greater than or equal to the observed
  return(r)
})

x_nonnectar_tb = x_nonnectar %>% as_tibble(rownames = "pathway_hsa") %>% 
  mutate(p_val = (value+1)/(N_random_samples+1)) %>%  #  (r+1)/(n+1)
  left_join(custom_anno_id) %>% 
  left_join(pathway_acc_con_sum2)
```

    ## Joining, by = "pathway_hsa"
    ## Joining, by = c("pathway_hsa", "description_hsa")

``` r
head(x_nonnectar_tb)
```

    ## # A tibble: 6 × 5
    ##   pathway_hsa   value  p_val description_hsa                          total_cn…¹
    ##   <chr>         <int>  <dbl> <chr>                                         <int>
    ## 1 path:hsa00010  2055 0.0206 Glycolysis / Gluconeogenesis                    689
    ## 2 path:hsa00020  4830 0.0483 Citrate cycle (TCA cycle)                       121
    ## 3 path:hsa00030  2850 0.0285 Pentose phosphate pathway                       459
    ## 4 path:hsa00040 29226 0.292  Pentose and glucuronate interconversions        846
    ## 5 path:hsa00051 18032 0.180  Fructose and mannose metabolism                 957
    ## 6 path:hsa00052 89574 0.896  Galactose metabolism                            777
    ## # … with abbreviated variable name ¹​total_cnees

``` r
hist(x_nonnectar_tb$p_val)
```

![](convergenceLv_v3_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
hist(x_nonnectar_tb$total_cnees)
```

![](convergenceLv_v3_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
ggplot(x_nonnectar_tb, aes(x = total_cnees, y = p_val)) +
  geom_point() 
```

![](convergenceLv_v3_files/figure-gfm/unnamed-chunk-7-3.png)<!-- -->

``` r
ggplot(x_nonnectar_tb, aes(x = total_cnees, y = p_val)) +
  geom_point() +
  scale_x_continuous(limits = c(NA, 2500))
```

    ## Warning: Removed 95 rows containing missing values (geom_point).

![](convergenceLv_v3_files/figure-gfm/unnamed-chunk-7-4.png)<!-- -->

``` r
x_nonnectar_tb_sig = x_nonnectar_tb %>% 
  filter(p_val<=0.05)
nrow(x_nonnectar_tb_sig) # p_val<= 0.05: 156; p_val<= 0.01: 96
```

    ## [1] 156

## comparison

``` r
ggvenn(list('nectar' = x_nectar_tb_sig$pathway_hsa,
            'non-nectar' = x_nonnectar_tb_sig$pathway_hsa))
```

![](convergenceLv_v3_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
onlyInNect = x_nectar_tb_sig %>% 
  filter(! pathway_hsa %in% x_nonnectar_tb_sig$pathway_hsa ) 
print(onlyInNect)
```

    ## # A tibble: 68 × 5
    ##    pathway_hsa   value     p_val description_hsa                         total…¹
    ##    <chr>         <int>     <dbl> <chr>                                     <int>
    ##  1 path:hsa00051   694 0.00695   Fructose and mannose metabolism             957
    ##  2 path:hsa00052   956 0.00957   Galactose metabolism                        777
    ##  3 path:hsa00071  3795 0.0380    Fatty acid degradation                      525
    ##  4 path:hsa00120  1312 0.0131    Primary bile acid biosynthesis              192
    ##  5 path:hsa00220   201 0.00202   Arginine biosynthesis                       254
    ##  6 path:hsa00260  2681 0.0268    Glycine, serine and threonine metaboli…     898
    ##  7 path:hsa00330     1 0.0000200 Arginine and proline metabolism             972
    ##  8 path:hsa00340   529 0.00530   Histidine metabolism                        623
    ##  9 path:hsa00380    38 0.000390  Tryptophan metabolism                      1075
    ## 10 path:hsa00400   725 0.00726   Phenylalanine, tyrosine and tryptophan…      43
    ## # … with 58 more rows, and abbreviated variable name ¹​total_cnees
