Methodological repetition of Zhu *et al.* study
================

## Introduction

In the study by [Zhu *et
al.*](https://gut.bmj.com/content/early/2023/01/23/gutjnl-2022-329307-0),
the authors present Mendelian Randomization (MR) analysis utilizing SNPs
of Heamorrhoidal disease (HEM) GWAS as exposure and irritable bowel
syndrome (IBS) as well as diverticular disease (DIV) GWAS as outcomes
and show that HEM significantly decreased the incidence of DIV as well
as of IBS. Here, we aim to methodologically repeat their analysis using
the same datasets, methods and parameters.

 

## Methods

### Instrument variable (IV) selection

Thresholds used by Zhu *et al.*:

1)  Association with exposure (HEM GWAS) - P \< 5e<sup>-10</sup>;
2)  Pruning for linkage disequilibrium (LD) - R<sup>2</sup> \< 0.001
    within 10kb (1000G EUR);
3)  No association with outcome (DIV and IBS GWAS) - P \< 0.05

 

### Mendelian Randomisation methods

MR methods used by Zhu *et al.*:

1)  Inverse variance weighted (IVW);
2)  MR Egger;
3)  Weighted median;
4)  Simple mode methods

 

### GWAS summary statistics data

Used in the Zhu *et al.* analysis:

- HEM ([Zheng *et al.*](https://pubmed.ncbi.nlm.nih.gov/33888516/)) -
  GCST90014033
- DIV ([Schafmayer *et
  al.*](https://pubmed.ncbi.nlm.nih.gov/30661054/)) - GCST008105
- DIV ([Dönertaş *et al.*](https://pubmed.ncbi.nlm.nih.gov/33959723/)) -
  GCST90038682
- IBS ([Dönertaş *et al.*](https://pubmed.ncbi.nlm.nih.gov/33959723/)) -
  GCST90038626
- IBS ([Jang et *al.*](https://pubmed.ncbi.nlm.nih.gov/34737426/)) -
  GCST90044173

 

### Packages

MR analyses were performed using R packages
[`TwoSampleMR`](https://mrcieu.github.io/TwoSampleMR/) and
[`tidyverse`](https://www.tidyverse.org).

``` r
library(tidyverse)
library(TwoSampleMR)
```

 

## Analysis and Results

 

### Data download

``` r
ebi_ftp <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"

# DIV - GCST008105
# https://pubmed.ncbi.nlm.nih.gov/30661054/
download.file( 
  url = paste0(ebi_ftp,"GCST008001-GCST009000/GCST008105/GWAS_summary_1-23.dosages.maf_0.01.info_0.4.txt.gz"),
  destfile = "../data/GWAS_summary_1-23.dosages.maf_0.01.info_0.4.txt.gz",
  method = "curl"
)

# DIV - GCST90038682
# https://pubmed.ncbi.nlm.nih.gov/33959723/
download.file( 
  url = paste0(ebi_ftp,"GCST90038001-GCST90039000/GCST90038682/GCST90038682_buildGRCh37.tsv"),
  destfile = "../data/GCST90038682_buildGRCh37.tsv",
  method = "curl"
)
# compress
R.utils::gzip(
  "../data/GCST90038682_buildGRCh37.tsv",
  overwrite=TRUE
)

# IBS - GCST90038626
# https://pubmed.ncbi.nlm.nih.gov/33959723/
download.file( 
  url = paste0(ebi_ftp,"GCST90038001-GCST90039000/GCST90038626/GCST90038626_buildGRCh37.tsv"),
  destfile = "../data/GCST90038626_buildGRCh37.tsv",
  method = "curl"
)
# compress
R.utils::gzip(
  "../data/GCST90038626_buildGRCh37.tsv",
  overwrite=TRUE
)

# IBS - GCST90044173
# https://pubmed.ncbi.nlm.nih.gov/34737426/
download.file( 
  url = paste0(ebi_ftp,"GCST90044001-GCST90045000/GCST90044173/GCST90044173_buildGRCh37.tsv.gz"),
  destfile = "../data/GCST90044173_buildGRCh37.tsv.gz",
  method = "curl"
)

# HEM - GCST90014033
# https://pubmed.ncbi.nlm.nih.gov/33888516/
# can be accessed directly from TwoSampleMR
```

 

### HEM (exposure) and DIV (outcome)

 

##### Obtain IVs from GWAS of exposure (HEM)

Here we obtain SNPs from HEM GWAS using clumping procedure with
thresholds:

``` r
# from OpenGWAS db
# thresholds as defined by Zhu et al.
HEM_exp_dat <- extract_instruments(
  "ebi-a-GCST90014033",
  p1 = 5e-10,
  r2 = 0.001,
  kb = 10
  ) %>% mutate(
    exposure = "HEM (Zheng et al.)"
  )

# number of SNPs after pruning
nrow(HEM_exp_dat)
```

    ## [1] 318

 

##### Obtain effects of IVs on outcome (DIV)

Here we import DIV GWAS datasets (Schafmayer *et al.* and Dönertaş *et
al.*) and filter the effects of SNPs, which will be used as IVs.

``` r
# Schafmayer et al. data
DIV_1_out_dat <- read_outcome_data(
    snps = HEM_exp_dat$SNP,
    filename = "../data/GWAS_summary_1-23.dosages.maf_0.01.info_0.4.txt.gz",
    sep = " ",
    snp_col = "SNP",
    chr_col = "CHR", 
    pos_col = "BP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "ALLELE1",
    other_allele_col = "ALLELE0",
    eaf_col = "A1FREQ",
    pval_col = "P"
    ) %>% 
  # add outcome name
  mutate(
    outcome = "DIV (Schafmayer et al.)"
  )

# Dönertaş et al. data
DIV_2_out_dat <- read_outcome_data(
    snps = HEM_exp_dat$SNP,
    filename = "../data/GCST90038682_buildGRCh37.tsv.gz",
    sep = "\t",
    snp_col = "variant_id",
    chr_col = "chromosome", 
    pos_col = "base_pair_location",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value"
    ) %>% 
  # add outcome name
  mutate(
    outcome = "DIV (Dönertaş et al.)"
  )

DIV_out_dat_list <- list(
  "DIV (Schafmayer et al.)" = DIV_1_out_dat,
  "DIV (Dönertaş et al.)" = DIV_2_out_dat
)

# check if any SNPs are missing
missing <- DIV_out_dat_list %>% 
  lapply(., function(x){
    x %>% filter(!SNP %in% HEM_exp_dat$SNP) %>% 
  pull(SNP)
  })
missing
```

    ## $`DIV (Schafmayer et al.)`
    ## character(0)
    ## 
    ## $`DIV (Dönertaş et al.)`
    ## character(0)

 

#### Harmonize exposure and outcome data

Here we harmonize SNPs, so the effects are relative to the same allele.
After harmonization, we then exclude SNPs, which have an association
with the outcome (P \< 0.05) as described by Zhu *et al.*

``` r
# harmonize data and exclude outcome-associated SNPs (p<0.05)
# as described by Zhu et al.
HEM_DIV_dat_list <- lapply(DIV_out_dat_list, function(x){
  harmonise_data(HEM_exp_dat, x) %>% 
    mutate(
      mr_keep = case_when(
        pval.outcome < 0.05 ~ FALSE,
        TRUE ~ mr_keep
        )
      )
  })

# number of SNPs retained for MR analysis
lapply(HEM_DIV_dat_list, function(x){
  x %>% filter(mr_keep == TRUE) %>% 
    nrow()
})
```

    ## $`DIV (Schafmayer et al.)`
    ## [1] 184
    ## 
    ## $`DIV (Dönertaş et al.)`
    ## [1] 255

 

#### Test for horizontal pleiotropy

Here we use MR Egger intercept test to evaluate if selected IVs are not
pleiotropic, thus have no effect on disease outside of its effect on the
exposure.

``` r
# test for horizontal pleiotropy
lapply(HEM_DIV_dat_list, mr_pleiotropy_test) %>% 
  bind_rows() %>% select(-1:-2)
```

    ##                   outcome           exposure egger_intercept           se
    ## 1 DIV (Schafmayer et al.) HEM (Zheng et al.)   -9.691874e-05 1.966784e-04
    ## 2   DIV (Dönertaş et al.) HEM (Zheng et al.)    5.443598e-05 4.952061e-05
    ##        pval
    ## 1 0.6227636
    ## 2 0.2727000

The results show that selected variants are unlikely to be horizontally
pleiotropic (P \> 0.05).

 

#### Perform MR analyses

Here, using the retained IVs, we perform MR by utilizing IVW, MR Egger,
weighted median and simple mode methods (the same as used by Zhu *et
al.*). In the analysis by Zhu et al, the authors call MR analysis using
DIV Schafmayer *et al.* and Dönertaş *et al.* datasets as “Set of
experiment” and “Set of verification”, respectively.

``` r
# perform MR
HEM_DIV_list_res <- lapply(
  HEM_DIV_dat_list,
  mr,
  method_list = list(
    "mr_ivw", "mr_egger_regression",
    "mr_weighted_median", "mr_simple_mode"
    )
  )

# calculate odds ratio
HEM_DIV_list_res <- lapply(HEM_DIV_list_res, function(x){
  generate_odds_ratios(x)
})

# results
HEM_DIV_list_res %>% bind_rows() %>% 
  select(-1:-2)
```

    ##                   outcome           exposure                    method nsnp
    ## 1 DIV (Schafmayer et al.) HEM (Zheng et al.) Inverse variance weighted  184
    ## 2 DIV (Schafmayer et al.) HEM (Zheng et al.)                  MR Egger  184
    ## 3 DIV (Schafmayer et al.) HEM (Zheng et al.)           Weighted median  184
    ## 4 DIV (Schafmayer et al.) HEM (Zheng et al.)               Simple mode  184
    ## 5   DIV (Dönertaş et al.) HEM (Zheng et al.) Inverse variance weighted  255
    ## 6   DIV (Dönertaş et al.) HEM (Zheng et al.)                  MR Egger  255
    ## 7   DIV (Dönertaş et al.) HEM (Zheng et al.)           Weighted median  255
    ## 8   DIV (Dönertaş et al.) HEM (Zheng et al.)               Simple mode  255
    ##               b           se         pval        lo_ci         up_ci        or
    ## 1 -0.0032271348 0.0014642373 2.752647e-02 -0.006097040 -0.0003572297 0.9967781
    ## 2 -0.0002652414 0.0061871070 9.658521e-01 -0.012391971  0.0118614883 0.9997348
    ## 3 -0.0047854467 0.0020659700 2.054083e-02 -0.008834748 -0.0007361456 0.9952260
    ## 4 -0.0121621658 0.0066899660 7.070379e-02 -0.025274499  0.0009501676 0.9879115
    ## 5 -0.0016096753 0.0004088385 8.244083e-05 -0.002410999 -0.0008083519 0.9983916
    ## 6 -0.0030362770 0.0013606594 2.652615e-02 -0.005703169 -0.0003693847 0.9969683
    ## 7 -0.0041770311 0.0006105863 7.863324e-12 -0.005373780 -0.0029802819 0.9958317
    ## 8 -0.0049921803 0.0017396259 4.454647e-03 -0.008401847 -0.0015825135 0.9950203
    ##    or_lci95  or_uci95
    ## 1 0.9939215 0.9996428
    ## 2 0.9876845 1.0119321
    ## 3 0.9912042 0.9992641
    ## 4 0.9750422 1.0009506
    ## 5 0.9975919 0.9991920
    ## 6 0.9943131 0.9996307
    ## 7 0.9946406 0.9970242
    ## 8 0.9916333 0.9984187

While results of IWV for the *set of experiment* and *set of
verification* were statistically significant (P = 0.0275 and P =
0.0000824, respectively) and the bxy values were negative as in the
study by Zhu et al., the obtained values were different than the ones in
the publication. Additionally, we were not able to recover exact number
of IVs, even after retaining all palindromic variants
\[`harmonise_data(..., action = 1)`\] during data harmonization step.

 

``` r
# plot results
HEM_DIV_list_plot <- lapply(seq(HEM_DIV_list_res), function(i){
  mr_scatter_plot(HEM_DIV_list_res[[i]], HEM_DIV_dat_list[[i]])
})
names(HEM_DIV_list_plot) <- names(HEM_DIV_list_res)

lapply(HEM_DIV_list_plot, function(x){
  x[[1]]
})
```

    ## $`DIV (Schafmayer et al.)`

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_repetition_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

    ## 
    ## $`DIV (Dönertaş et al.)`

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_repetition_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

 

### HEM (exposure) and IBS (outcome)

 

##### Obtain effects of IVs on outcome (IBS)

Here we import IBS GWAS datasets (Dönertaş *et al.* and Jiang *et al.*)
and filter the effects of SNPs, which will be used as IVs.

``` r
# list files to IBS studies
IBS_file_list <- list(
  "IBS (Dönertaş et al.)" = "../data/GCST90038626_buildGRCh37.tsv.gz",
  "IBS (Jiang et al.)" = "../data/GCST90044173_buildGRCh37.tsv.gz"
)

# get effects of instruments on outcome
IBS_out_dat_list <- lapply(IBS_file_list, function(x){
  read_outcome_data(
    snps = HEM_exp_dat$SNP,
    filename = x,
    sep = "\t",
    snp_col = "variant_id",
    chr_col = "chromosome", 
    pos_col = "base_pair_location",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value"
    )
  })

# add outcome name
for(i in seq(IBS_out_dat_list)){
  IBS_out_dat_list[[i]] <- IBS_out_dat_list[[i]] %>% 
    mutate(
      outcome = names(IBS_out_dat_list[i])
    )
}

# check if any SNPs are missing
missing <- IBS_out_dat_list %>% 
  lapply(., function(x){
    x %>% filter(!SNP %in% HEM_exp_dat$SNP) %>% 
  pull(SNP)
  })
missing
```

    ## $`IBS (Dönertaş et al.)`
    ## character(0)
    ## 
    ## $`IBS (Jiang et al.)`
    ## character(0)

 

#### Harmonize exposure and outcome data

As previously, we harmonize and then exclude SNPs, which have an
association with the outcome (P \< 0.05) as described by Zhu *et al.*

``` r
# harmonize data and exclude outcome-associated SNPs (p<0.05)
# as described by Zhu et al.
HEM_IBS_dat_list <- lapply(IBS_out_dat_list, function(x){
  harmonise_data(HEM_exp_dat, x) %>% 
    mutate(
      mr_keep = case_when(
        pval.outcome < 0.05 ~ FALSE,
        TRUE ~ mr_keep
        )
      )
  })

# number of SNPs retained for MR analysis
lapply(HEM_IBS_dat_list, function(x){
  x %>% filter(mr_keep == TRUE) %>% 
    nrow()
})
```

    ## $`IBS (Dönertaş et al.)`
    ## [1] 278
    ## 
    ## $`IBS (Jiang et al.)`
    ## [1] 290

 

#### Test for horizontal pleiotropy

As before, we use MR Egger intercept test to evaluate if selected IVs
are not horizontally pleiotropic.

``` r
# test for horizontal pleiotropy
lapply(HEM_IBS_dat_list, mr_pleiotropy_test) %>% 
  bind_rows() %>% select(-1:-2)
```

    ##                 outcome           exposure egger_intercept           se
    ## 1 IBS (Dönertaş et al.) HEM (Zheng et al.)    0.0002548583 6.400372e-05
    ## 2    IBS (Jiang et al.) HEM (Zheng et al.)    0.0389544589 7.638534e-03
    ##           pval
    ## 1 8.748380e-05
    ## 2 6.176597e-07

This time, the results show that the selected variants are likely to be
horizontally pleiotropic (P \< 0.05).

 

#### Perform MR analyses

As previously, using the retained IVs, we perform MR by utilizing IVW,
MR Egger, weighted median and simple mode methods (the same as used by
Zhu *et al.*). As in case for DIV, the authors of the Zhu *et al.* call
MR analysis using IBS Dönertaş *et al.* and Jiang *et al.* datasets as
“Set of experiment” and “Set of verification”, respectively.

``` r
# perform MR
HEM_IBS_list_res <- lapply(
  HEM_IBS_dat_list,
  mr,
  method_list = list(
    "mr_ivw", "mr_egger_regression",
    "mr_weighted_median", "mr_simple_mode"
    )
  )

# calculate odds ratio
HEM_IBS_list_res <- lapply(HEM_IBS_list_res, function(x){
  generate_odds_ratios(x)
})

# results
HEM_IBS_list_res %>% bind_rows() %>% 
  select(-1:-2)
```

    ##                 outcome           exposure                    method nsnp
    ## 1 IBS (Dönertaş et al.) HEM (Zheng et al.) Inverse variance weighted  278
    ## 2 IBS (Dönertaş et al.) HEM (Zheng et al.)                  MR Egger  278
    ## 3 IBS (Dönertaş et al.) HEM (Zheng et al.)           Weighted median  278
    ## 4 IBS (Dönertaş et al.) HEM (Zheng et al.)               Simple mode  278
    ## 5    IBS (Jiang et al.) HEM (Zheng et al.) Inverse variance weighted  290
    ## 6    IBS (Jiang et al.) HEM (Zheng et al.)                  MR Egger  290
    ## 7    IBS (Jiang et al.) HEM (Zheng et al.)           Weighted median  290
    ## 8    IBS (Jiang et al.) HEM (Zheng et al.)               Simple mode  290
    ##              b           se         pval        lo_ci         up_ci        or
    ## 1 -0.001891724 0.0005792522 1.091536e-03 -0.003027059 -0.0007563901 0.9981101
    ## 2 -0.008460040 0.0017454651 2.095793e-06 -0.011881152 -0.0050389287 0.9915756
    ## 3 -0.004727770 0.0008855887 9.368689e-08 -0.006463523 -0.0029920158 0.9952834
    ## 4 -0.006984177 0.0030819910 2.421539e-02 -0.013024880 -0.0009434751 0.9930402
    ## 5 -0.082729407 0.0673755047 2.194899e-01 -0.214785396  0.0493265819 0.9206002
    ## 6 -1.088852423 0.2083032936 3.310516e-07 -1.497126879 -0.6805779679 0.3366025
    ## 7 -0.449845216 0.0990530469 5.586661e-06 -0.643989188 -0.2557012441 0.6377269
    ## 8 -0.742966359 0.3146923310 1.889306e-02 -1.359763328 -0.1261693904 0.4757007
    ##    or_lci95  or_uci95
    ## 1 0.9969775 0.9992439
    ## 2 0.9881892 0.9949737
    ## 3 0.9935573 0.9970125
    ## 4 0.9870596 0.9990570
    ## 5 0.8067145 1.0505634
    ## 6 0.2237722 0.5063243
    ## 7 0.5251931 0.7743733
    ## 8 0.2567215 0.8814655

Similarity to DIV, the results of IWV for the *set of experiment* and
*set of verification* were different from the ones published by Zhu *et
al.* As in the case for DIV, we were not able to recover exact number of
IVs, even after retaining all palindromic variants
\[`harmonise_data(..., action = 1)`\] during data harmonization step.

 

``` r
# plot results
HEM_IBS_list_plot <- lapply(seq(HEM_IBS_list_res), function(i){
  mr_scatter_plot(HEM_IBS_list_res[[i]], HEM_IBS_dat_list[[i]])
})
names(HEM_IBS_list_plot) <- names(HEM_IBS_list_res)

lapply(HEM_IBS_list_plot, function(x){
  x[[1]]
})
```

    ## $`IBS (Dönertaş et al.)`

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_repetition_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

    ## 
    ## $`IBS (Jiang et al.)`

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_repetition_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

 

In conclusion, by following the instructions provided by Zhu *et al.*,
we could not reproduce the exact results as described in their letter.

 

``` r
sessionInfo()
```

    ## R version 4.2.1 (2022-06-23)
    ## Platform: aarch64-apple-darwin20 (64-bit)
    ## Running under: macOS Monterey 12.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] TwoSampleMR_0.5.6 forcats_0.5.2     stringr_1.4.1     dplyr_1.1.2      
    ##  [5] purrr_0.3.4       readr_2.1.2       tidyr_1.2.1       tibble_3.2.1     
    ##  [9] ggplot2_3.4.2     tidyverse_1.3.2  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] httr_1.4.4          jsonlite_1.8.0      splines_4.2.1      
    ##  [4] foreach_1.5.2       R.utils_2.12.2      modelr_0.1.9       
    ##  [7] assertthat_0.2.1    highr_0.9           googlesheets4_1.0.1
    ## [10] cellranger_1.1.0    yaml_2.3.5          pillar_1.9.0       
    ## [13] backports_1.4.1     lattice_0.20-45     glue_1.6.2         
    ## [16] digest_0.6.29       rvest_1.0.3         colorspace_2.0-3   
    ## [19] R.oo_1.25.0         htmltools_0.5.5     Matrix_1.5-1       
    ## [22] plyr_1.8.7          pkgconfig_2.0.3     broom_1.0.1        
    ## [25] haven_2.5.1         scales_1.2.1        tzdb_0.3.0         
    ## [28] timechange_0.2.0    googledrive_2.0.0   farver_2.1.1       
    ## [31] generics_0.1.3      ellipsis_0.3.2      withr_2.5.0        
    ## [34] cli_3.4.1           survival_3.4-0      magrittr_2.0.3     
    ## [37] crayon_1.5.1        readxl_1.4.1        ieugwasr_0.1.5     
    ## [40] evaluate_0.16       R.methodsS3_1.8.2   fs_1.5.2           
    ## [43] fansi_1.0.3         xml2_1.3.3          tools_4.2.1        
    ## [46] data.table_1.14.2   hms_1.1.2           gargle_1.2.1       
    ## [49] lifecycle_1.0.3     munsell_0.5.0       reprex_2.0.2       
    ## [52] glmnet_4.1-4        mr.raps_0.2         compiler_4.2.1     
    ## [55] rlang_1.1.1         grid_4.2.1          iterators_1.0.14   
    ## [58] rstudioapi_0.14     labeling_0.4.2      rmarkdown_2.16     
    ## [61] gtable_0.3.1        codetools_0.2-18    DBI_1.1.3          
    ## [64] curl_4.3.2          R6_2.5.1            lubridate_1.9.2    
    ## [67] knitr_1.40          fastmap_1.1.0       utf8_1.2.2         
    ## [70] nortest_1.0-4       shape_1.4.6         stringi_1.7.8      
    ## [73] Rcpp_1.0.9          vctrs_0.6.3         dbplyr_2.2.1       
    ## [76] tidyselect_1.2.0    xfun_0.33
