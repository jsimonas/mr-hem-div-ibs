Methodological refinement of Zhu *et al.* study
================

## Introduction

In the study by [Zhu *et
al.*](https://gut.bmj.com/content/early/2023/01/23/gutjnl-2022-329307-0),
the authors used correlated (in linkage disequilibrium) genetic variants
as IVs, which violates assumptions of the IWV method [(Burgess and
Thompson,
2017)](https://link.springer.com/article/10.1007/s10654-017-0255-x).
Additionally, to remove horizontal pleiotropy, the authors used a
threshold of P\<0.05 to test for the association to the outcome on a
genome-wide scale. Due to the high false positive rate of the P\<0.05,
this may lead to the removal of vertically pleiotropic variants.

Here, we refined MR analyses of Zhu *et al.* study, by employing a
standard workflow for IV selection to ensure that the variants are
uncorrelated as well as not horizontally pleiotropic.

 

## Methods

### Instrument variable (IV) selection

Refined thresholds used:

1)  Association with exposure (HEM GWAS) - P \< 5e<sup>-8</sup>;
2)  Clumping for linkage disequilibrium (LD) - r<sup>2</sup> \< 0.001
    within 1 Mb distance (1000G EUR);
3)  MR-PRESSO outlier test (P\<0.05)

 

### Mendelian Randomisation methods

MR methods used:

1)  Inverse variance weighted (IVW);
2)  MR Egger;
3)  Weighted median;
4)  Weighted mode;
5)  Simple mode

 

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
  destfile = "data/GCST90038682_buildGRCh37.tsv",
  method = "curl"
)
# compress
R.utils::gzip(
  "data/GCST90038682_buildGRCh37.tsv",
  overwrite=TRUE
)

# IBS - GCST90038626
# https://pubmed.ncbi.nlm.nih.gov/33959723/
download.file( 
  url = paste0(ebi_ftp,"GCST90038001-GCST90039000/GCST90038626/GCST90038626_buildGRCh37.tsv"),
  destfile = "data/GCST90038626_buildGRCh37.tsv",
  method = "curl"
)
# compress
R.utils::gzip(
  "data/GCST90038626_buildGRCh37.tsv",
  overwrite=TRUE
)

# IBS - GCST90044173
# https://pubmed.ncbi.nlm.nih.gov/34737426/
download.file( 
  url = paste0(ebi_ftp,"GCST90044001-GCST90045000/GCST90044173/GCST90044173_buildGRCh37.tsv.gz"),
  destfile = "data/GCST90044173_buildGRCh37.tsv.gz",
  method = "curl"
)

# HEM - GCST90014033
# https://pubmed.ncbi.nlm.nih.gov/33888516/
# can be accessed directly from TwoSampleMR
```

 

### HEM (exposure) and DIV (outcome)

 

##### Obtain IVs from GWAS of exposure (HEM)

Here we obtain SNPs from HEM GWAS using clumping procedure with refined
thresholds:

``` r
# from OpenGWAS db
# thresholds as defined by Zhu et al.
HEM_exp_dat <- extract_instruments(
  "ebi-a-GCST90014033",
  p1 = 5e-8,
  r2 = 0.001,
  kb = 1000
  ) %>% mutate(
    exposure = "HEM (Zheng et al.)"
  )

# number of SNPs after clumping
nrow(HEM_exp_dat)
```

    ## [1] 101

 

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

``` r
# harmonize data
HEM_DIV_dat_list <- lapply(DIV_out_dat_list, function(x){
  harmonise_data(HEM_exp_dat, x)
  })

# number of SNPs retained for MR analysis
lapply(HEM_DIV_dat_list, function(x){
  x %>% filter(mr_keep == TRUE) %>% 
    nrow()
})
```

    ## $`DIV (Schafmayer et al.)`
    ## [1] 96
    ## 
    ## $`DIV (Dönertaş et al.)`
    ## [1] 96

 

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
    ## 1 DIV (Schafmayer et al.) HEM (Zheng et al.)   -0.0005415101 0.0004250404
    ## 2   DIV (Dönertaş et al.) HEM (Zheng et al.)   -0.0001468679 0.0001109848
    ##        pval
    ## 1 0.2057971
    ## 2 0.1889405

The results show that selected variants are no likely to be horizontally
pleiotropic (P \< 0.05).

 

Here we use MR-PRESSO to detect and remove ouliers and then again run MR
Egger intercept test to evaluate horizontal pleiotropy.

``` r
# set random seed
set.seed(0)

# run mr-presso to perform correction of
# horizontal pleiotropy via outlier removal
HEM_DIV_presso_list <- lapply(
  HEM_DIV_dat_list,
  run_mr_presso,
  NbDistribution = 2000,
  SignifThreshold = 0.1
  )

# global test results
lapply(HEM_DIV_presso_list, function(x){
  x[[1]]$`MR-PRESSO results`$`Global Test`
})
```

    ## $`DIV (Schafmayer et al.)`
    ## $`DIV (Schafmayer et al.)`$RSSobs
    ## [1] 431.7767
    ## 
    ## $`DIV (Schafmayer et al.)`$Pvalue
    ## [1] "<5e-04"
    ## 
    ## 
    ## $`DIV (Dönertaş et al.)`
    ## $`DIV (Dönertaş et al.)`$RSSobs
    ## [1] 189.4332
    ## 
    ## $`DIV (Dönertaş et al.)`$Pvalue
    ## [1] "<5e-04"

``` r
# HP outliers
lapply(names(HEM_DIV_dat_list), function(i){
  HEM_DIV_dat_list[[i]] %>% 
    filter(
      row_number() %in% HEM_DIV_presso_list[[i]][[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
      ) %>% 
  select(1:9)
})
```

    ## [[1]]
    ##           SNP effect_allele.exposure other_allele.exposure
    ## 1  rs11942410                      T                     C
    ## 2  rs13017210                      T                     A
    ## 3     rs13632                      G                     A
    ## 4  rs17293632                      T                     C
    ## 5   rs2180811                      A                     T
    ## 6      rs3253                      T                     C
    ## 7   rs3757582                      C                     T
    ## 8   rs4556017                      T                     C
    ## 9  rs62368263                      C                     T
    ## 10  rs7423637                      T                     A
    ## 11  rs8106090                      G                     A
    ##    effect_allele.outcome other_allele.outcome beta.exposure beta.outcome
    ## 1                      T                    C       -0.0267  0.002843160
    ## 2                      T                    A       -0.0259  0.001617480
    ## 3                      G                    A       -0.0304  0.000708936
    ## 4                      T                    C        0.0543  0.001671520
    ## 5                      A                    T        0.0278  0.000941901
    ## 6                      T                    C        0.0291  0.001643780
    ## 7                      C                    T        0.0558  0.007060800
    ## 8                      T                    C       -0.0569  0.002780150
    ## 9                      C                    T       -0.0441  0.000926076
    ## 10                     T                    A        0.0244  0.001650880
    ## 11                     G                    A        0.0217  0.000270501
    ##    eaf.exposure eaf.outcome
    ## 1        0.2552    0.260126
    ## 2        0.3859    0.401810
    ## 3        0.7639    0.768868
    ## 4        0.2373    0.236850
    ## 5        0.4831    0.477374
    ## 6        0.3202    0.313720
    ## 7        0.0560    0.061507
    ## 8        0.8556    0.852616
    ## 9        0.1426    0.132724
    ## 10       0.3086    0.306306
    ## 11       0.4920    0.485250
    ## 
    ## [[2]]
    ##          SNP effect_allele.exposure other_allele.exposure effect_allele.outcome
    ## 1 rs11942410                      T                     C                     T
    ## 2  rs2180811                      A                     T                     A
    ## 3     rs3253                      T                     C                     T
    ## 4  rs7423637                      T                     A                     T
    ##   other_allele.outcome beta.exposure beta.outcome eaf.exposure eaf.outcome
    ## 1                    C       -0.0267  7.96604e-04       0.2552    0.259579
    ## 2                    T        0.0278  6.12777e-05       0.4831    0.488814
    ## 3                    C        0.0291  2.55194e-04       0.3202    0.319748
    ## 4                    A        0.0244  4.90248e-04       0.3086    0.303789

``` r
# remove HP outliers
HEM_DIV_dat_list_adj <- lapply(names(HEM_DIV_dat_list), function(i){
    HEM_DIV_dat_list[[i]] %>% 
    filter(
      ! row_number() %in% HEM_DIV_presso_list[[i]][[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
      ) 
})
names(HEM_DIV_dat_list_adj) <- names(HEM_DIV_dat_list)

# re-test for horizontal pleiotropy
lapply(HEM_DIV_dat_list_adj, mr_pleiotropy_test) %>% 
  bind_rows() %>% select(-1:-2)
```

    ##                   outcome           exposure egger_intercept           se
    ## 1 DIV (Schafmayer et al.) HEM (Zheng et al.)   -0.0005049278 0.0004150861
    ## 2   DIV (Dönertaş et al.) HEM (Zheng et al.)   -0.0001507380 0.0001084789
    ##        pval
    ## 1 0.2272240
    ## 2 0.1680529

The results show that selected variants are not likely to be
horizontally pleiotropic (P \> 0.05).

 

#### Perform MR analyses

Here, using the retained IVs, we perform MR by utilizing IVW, MR Egger,
weighted median, simple mode and Weighted mode methods (default of the
`TwoSampleMR` `mr()` function). In the analysis by Zhu et al, the
authors call MR analysis using DIV Schafmayer *et al.* and Dönertaş *et
al.* datasets as “Set of experiment” and “Set of verification”,
respectively.

``` r
# perform MR
HEM_DIV_list_res <- lapply(
  HEM_DIV_dat_list_adj, mr
  )

# calculate odds ratio
HEM_DIV_list_res <- lapply(HEM_DIV_list_res, function(x){
  generate_odds_ratios(x)
})

# results
HEM_DIV_list_res %>% bind_rows() %>% 
  select(-1:-2)
```

    ##                    outcome           exposure                    method nsnp
    ## 1  DIV (Schafmayer et al.) HEM (Zheng et al.)                  MR Egger   86
    ## 2  DIV (Schafmayer et al.) HEM (Zheng et al.)           Weighted median   86
    ## 3  DIV (Schafmayer et al.) HEM (Zheng et al.) Inverse variance weighted   86
    ## 4  DIV (Schafmayer et al.) HEM (Zheng et al.)               Simple mode   86
    ## 5  DIV (Schafmayer et al.) HEM (Zheng et al.)             Weighted mode   86
    ## 6    DIV (Dönertaş et al.) HEM (Zheng et al.)                  MR Egger   93
    ## 7    DIV (Dönertaş et al.) HEM (Zheng et al.)           Weighted median   93
    ## 8    DIV (Dönertaş et al.) HEM (Zheng et al.) Inverse variance weighted   93
    ## 9    DIV (Dönertaş et al.) HEM (Zheng et al.)               Simple mode   93
    ## 10   DIV (Dönertaş et al.) HEM (Zheng et al.)             Weighted mode   93
    ##                b          se       pval         lo_ci       up_ci        or
    ## 1   0.0217030325 0.012788052 0.09337390 -0.0033615493 0.046767614 1.0219403
    ## 2   0.0014258498 0.003387447 0.67381226 -0.0052135459 0.008065246 1.0014269
    ## 3   0.0068801312 0.003890200 0.07696392 -0.0007446603 0.014504923 1.0069039
    ## 4   0.0013213924 0.007779325 0.86552425 -0.0139260849 0.016568870 1.0013223
    ## 5  -0.0006624167 0.006774832 0.92234028 -0.0139410879 0.012616255 0.9993378
    ## 6   0.0051723050 0.003294491 0.11988913 -0.0012848974 0.011629507 1.0051857
    ## 7  -0.0001482294 0.001238514 0.90473405 -0.0025757166 0.002279258 0.9998518
    ## 8   0.0008207313 0.001028214 0.42474845 -0.0011945688 0.002836032 1.0008211
    ## 9   0.0041128359 0.003287021 0.21402120 -0.0023297262 0.010555398 1.0041213
    ## 10 -0.0036857831 0.003206311 0.25331297 -0.0099701529 0.002598587 0.9963210
    ##     or_lci95 or_uci95
    ## 1  0.9966441 1.047878
    ## 2  0.9948000 1.008098
    ## 3  0.9992556 1.014611
    ## 4  0.9861704 1.016707
    ## 5  0.9861556 1.012696
    ## 6  0.9987159 1.011697
    ## 7  0.9974276 1.002282
    ## 8  0.9988061 1.002840
    ## 9  0.9976730 1.010611
    ## 10 0.9900794 1.002602

``` r
# plot results
HEM_DIV_list_plot <- lapply(seq(HEM_DIV_list_res), function(i){
  mr_scatter_plot(HEM_DIV_list_res[[i]], HEM_DIV_dat_list_adj[[i]])
})
names(HEM_DIV_list_plot) <- names(HEM_DIV_list_res)

lapply(HEM_DIV_list_plot, function(x){
  x[[1]]
})
```

    ## $`DIV (Schafmayer et al.)`

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_refinement_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

    ## 
    ## $`DIV (Dönertaş et al.)`

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_refinement_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

No statistically significant associations were identified in both
analyses.

 

#### Perform sensitivity analyses

Here, we perform leave-one-out analysis to test if association was not
driven by a single variant.

``` r
# leave one out analysis
HEM_DIV_list_loo <- lapply(
  HEM_DIV_dat_list_adj, mr_leaveoneout
  )

# plot results
lapply(HEM_DIV_list_loo, function(x){
  p <- mr_leaveoneout_plot(x)
  p[[1]] + coord_flip() +
  theme(
    axis.text.x = element_text(
      angle = 90, vjust=0.5, hjust = 1, size = 6
    ))
  })
```

    ## $`DIV (Schafmayer et al.)`

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_refinement_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

    ## 
    ## $`DIV (Dönertaş et al.)`

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_refinement_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->
 

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

As previously, we harmonize SNPs in exposure and outcome datasets.

``` r
# harmonize data
HEM_IBS_dat_list <- lapply(IBS_out_dat_list, function(x){
  harmonise_data(HEM_exp_dat, x)
  })

# number of SNPs retained for MR analysis
lapply(HEM_IBS_dat_list, function(x){
  x %>% filter(mr_keep == TRUE) %>% 
    nrow()
})
```

    ## $`IBS (Dönertaş et al.)`
    ## [1] 96
    ## 
    ## $`IBS (Jiang et al.)`
    ## [1] 96

 

#### Test for horizontal pleiotropy

As before, we use MR Egger intercept test to evaluate if selected IVs
are not horizontally pleiotropic.

``` r
# test for horizontal pleiotropy
lapply(HEM_IBS_dat_list, mr_pleiotropy_test) %>% 
  bind_rows() %>% select(-1:-2)
```

    ##                 outcome           exposure egger_intercept           se
    ## 1 IBS (Dönertaş et al.) HEM (Zheng et al.)    0.0001466182 0.0001289706
    ## 2    IBS (Jiang et al.) HEM (Zheng et al.)    0.0218143172 0.0167121964
    ##        pval
    ## 1 0.2584972
    ## 2 0.1949790

 

Here we use MR-PRESSO to detect and remove ouliers and then again run MR
Egger intercept test to evaluate horizontal pleiotropy.

``` r
# set random seed
set.seed(0)

# run mr-presso to perform correction of
# horizontal pleiotropy via outlier removal
HEM_IBS_presso_list <- lapply(
  HEM_IBS_dat_list,
  run_mr_presso,
  NbDistribution = 2000,
  SignifThreshold = 0.1
  )

# global test results
lapply(HEM_IBS_presso_list, function(x){
  x[[1]]$`MR-PRESSO results`$`Global Test`
})
```

    ## $`IBS (Dönertaş et al.)`
    ## $`IBS (Dönertaş et al.)`$RSSobs
    ## [1] 120.3473
    ## 
    ## $`IBS (Dönertaş et al.)`$Pvalue
    ## [1] 0.066
    ## 
    ## 
    ## $`IBS (Jiang et al.)`
    ## $`IBS (Jiang et al.)`$RSSobs
    ## [1] 139.6664
    ## 
    ## $`IBS (Jiang et al.)`$Pvalue
    ## [1] 0.0025

Both methods show that the selected variants are not likely to be
horizontally pleiotropic (P \< 0.05).

 

#### Perform MR analyses

As previously, using the retained IVs, we perform MR by utilizing
default methods of the `TwoSampleMR::mr()` function. As in case for DIV,
the authors of the Zhu *et al.* call MR analysis using IBS Dönertaş *et
al.* and Jiang *et al.* datasets as “Set of experiment” and “Set of
verification”, respectively.

``` r
# perform MR
HEM_IBS_list_res <- lapply(
  HEM_IBS_dat_list, mr
  )

# calculate odds ratio
HEM_IBS_list_res <- lapply(HEM_IBS_list_res, function(x){
  generate_odds_ratios(x)
})

# results
HEM_IBS_list_res %>% bind_rows() %>% 
  select(-1:-2)
```

    ##                  outcome           exposure                    method nsnp
    ## 1  IBS (Dönertaş et al.) HEM (Zheng et al.)                  MR Egger   96
    ## 2  IBS (Dönertaş et al.) HEM (Zheng et al.)           Weighted median   96
    ## 3  IBS (Dönertaş et al.) HEM (Zheng et al.) Inverse variance weighted   96
    ## 4  IBS (Dönertaş et al.) HEM (Zheng et al.)               Simple mode   96
    ## 5  IBS (Dönertaş et al.) HEM (Zheng et al.)             Weighted mode   96
    ## 6     IBS (Jiang et al.) HEM (Zheng et al.)                  MR Egger   96
    ## 7     IBS (Jiang et al.) HEM (Zheng et al.)           Weighted median   96
    ## 8     IBS (Jiang et al.) HEM (Zheng et al.) Inverse variance weighted   96
    ## 9     IBS (Jiang et al.) HEM (Zheng et al.)               Simple mode   96
    ## 10    IBS (Jiang et al.) HEM (Zheng et al.)             Weighted mode   96
    ##               b          se      pval         lo_ci       up_ci        or
    ## 1  -0.002368259 0.003939219 0.5491537 -0.0100891290 0.005352611 0.9976345
    ## 2   0.000740980 0.001677650 0.6587227 -0.0025472136 0.004029174 1.0007413
    ## 3   0.001891967 0.001215990 0.1197304 -0.0004913742 0.004275309 1.0018938
    ## 4  -0.001943580 0.003989713 0.6272758 -0.0097634171 0.005876257 0.9980583
    ## 5  -0.002145878 0.003501644 0.5414595 -0.0090091008 0.004717344 0.9978564
    ## 6  -0.715765306 0.510543791 0.1642188 -1.7164311358 0.284900524 0.4888179
    ## 7  -0.040396762 0.214121163 0.8503569 -0.4600742404 0.379280717 0.9604083
    ## 8  -0.081770340 0.157871183 0.6044889 -0.3911978595 0.227657180 0.9214836
    ## 9  -0.118585070 0.502639600 0.8139990 -1.1037586848 0.866588545 0.8881763
    ## 10 -0.175153941 0.450757544 0.6984587 -1.0586387276 0.708330845 0.8393278
    ##     or_lci95 or_uci95
    ## 1  0.9899616 1.005367
    ## 2  0.9974560 1.004037
    ## 3  0.9995087 1.004284
    ## 4  0.9902841 1.005894
    ## 5  0.9910314 1.004728
    ## 6  0.1797064 1.329630
    ## 7  0.6312368 1.461233
    ## 8  0.6762463 1.255655
    ## 9  0.3316223 2.378782
    ## 10 0.3469278 2.030599

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

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_refinement_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

    ## 
    ## $`IBS (Jiang et al.)`

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_refinement_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

As in case to DIV, no statistically significant associations were
identified in both datasets.

 

#### Perform sensitivity analyses

Again, we perform leave-one-out analysis to test if association was not
driven by a single variant.

``` r
# leave one out analysis
HEM_IBS_list_loo <- lapply(
  HEM_IBS_dat_list, mr_leaveoneout
  )

# plot results
lapply(HEM_IBS_list_loo, function(x){
  p <- mr_leaveoneout_plot(x)
  p[[1]] + coord_flip() +
  theme(
    axis.text.x = element_text(
      angle = 90, vjust=0.5, hjust = 1, size = 6
    ))
  })
```

    ## $`IBS (Dönertaş et al.)`

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_refinement_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

    ## 
    ## $`IBS (Jiang et al.)`

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_refinement_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->
 

### DIV (exposure) and HEM (outcome)

 

##### Obtain effects of IVs on outcome (HEM)

Here we import DIV GWAS datasets (Schafmayer *et al.* and Dönertaş *et
al.*) and use as exposure data.

``` r
# get effects of instruments on outcome
DIV_1_exp_dat <- read_exposure_data(
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
  filter(
    pval.exposure < 5e-08
  ) %>% 
  mutate(
    exposure = "DIV (Schafmayer et al.)"
  )

DIV_2_exp_dat <- read_exposure_data(
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
  filter(
    pval.exposure < 1e-05
  ) %>% 
  mutate(
    exposure = "DIV (Dönertaş et al.)"
  )

DIV_exp_dat_list <- list(
  "DIV (Schafmayer et al.)" = DIV_1_exp_dat,
  "DIV (Dönertaş et al.)" = DIV_2_exp_dat
)

# add exposure name
for(i in seq(DIV_exp_dat_list)){
  DIV_exp_dat_list[[i]] <- DIV_exp_dat_list[[i]] %>% 
    mutate(
      exposure = names(DIV_exp_dat_list[i])
    )
}

# clump instruments
DIV_exp_dat_list <- lapply(
  DIV_exp_dat_list,
  clump_data,
  clump_r2 = 0.001,
  clump_kb = 1000
  )

# number of SNPs
lapply(DIV_exp_dat_list, nrow)
```

    ## $`DIV (Schafmayer et al.)`
    ## [1] 53
    ## 
    ## $`DIV (Dönertaş et al.)`
    ## [1] 50

``` r
# get effects of instruments on outcome
HEM_out_dat_list <- lapply(DIV_exp_dat_list, function(x){
  extract_outcome_data(
    snps = x$SNP,
    outcomes = "ebi-a-GCST90014033"
    ) %>% 
    mutate(
      outcome = "HEM (Zheng et al.)"
      )
  })

# check if any SNPs are missing
missing <- lapply(seq(HEM_out_dat_list), function(i){
    DIV_exp_dat_list[[i]] %>% filter(!SNP %in% HEM_out_dat_list[[i]]$SNP) %>% 
    pull(SNP)
  })
names(missing) <- names(HEM_out_dat_list)
missing
```

    ## $`DIV (Schafmayer et al.)`
    ## [1] "rs7990"     "rs60869342" "rs3752946"  "rs77161253"
    ## 
    ## $`DIV (Dönertaş et al.)`
    ## [1] "rs13020164"  "rs185433580" "rs190098160" "rs56308199"  "rs11083450"

 

#### Harmonize exposure and outcome data

As previously, we harmonize SNPs in exposure and outcome datasets.

``` r
# harmonize the exposure and outcome data
DIV_HEM_dat_list <- lapply(seq(DIV_exp_dat_list), function(i){
  harmonise_data(DIV_exp_dat_list[[i]], HEM_out_dat_list[[i]])
})
names(DIV_HEM_dat_list) <- names(DIV_out_dat_list)

# number of SNPs retained for MR analysis
lapply(DIV_HEM_dat_list, function(x){
  x %>% filter(mr_keep == TRUE) %>% 
    nrow()
})
```

    ## $`DIV (Schafmayer et al.)`
    ## [1] 48
    ## 
    ## $`DIV (Dönertaş et al.)`
    ## [1] 45

 

#### Test for horizontal pleiotropy

As before, we use MR Egger intercept test to evaluate if selected IVs
are not horizontally pleiotropic.

``` r
# test for horizontal pleiotropy
lapply(DIV_HEM_dat_list, mr_pleiotropy_test) %>% 
  bind_rows() %>% select(-1:-2)
```

    ##              outcome                exposure egger_intercept          se
    ## 1 HEM (Zheng et al.) DIV (Schafmayer et al.)    -0.006833984 0.005715271
    ## 2 HEM (Zheng et al.)   DIV (Dönertaş et al.)     0.002582882 0.003428257
    ##        pval
    ## 1 0.2379254
    ## 2 0.4553093

 

Here we use MR-PRESSO to detect and remove outliers and then again run
MR Egger intercept test to evaluate horizontal pleiotropy.

``` r
# set random seed
set.seed(0)

# run mr-presso to perform correction of
# horizontal pleiotropy via outlier removal
DIV_HEM_presso_list <- lapply(
  DIV_HEM_dat_list,
  run_mr_presso,
  NbDistribution = 2000,
  SignifThreshold = 0.1
  )

# global test results
lapply(DIV_HEM_presso_list, function(x){
  x[[1]]$`MR-PRESSO results`$`Global Test`
})
```

    ## $`DIV (Schafmayer et al.)`
    ## $`DIV (Schafmayer et al.)`$RSSobs
    ## [1] 348.4132
    ## 
    ## $`DIV (Schafmayer et al.)`$Pvalue
    ## [1] "<5e-04"
    ## 
    ## 
    ## $`DIV (Dönertaş et al.)`
    ## $`DIV (Dönertaş et al.)`$RSSobs
    ## [1] 83.04019
    ## 
    ## $`DIV (Dönertaş et al.)`$Pvalue
    ## [1] 0.002

``` r
# HP outliers
lapply(names(DIV_HEM_dat_list), function(i){
  DIV_HEM_dat_list[[i]] %>% 
    filter(
      row_number() %in% DIV_HEM_presso_list[[i]][[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
      ) %>% 
  select(1:9)
})
```

    ## [[1]]
    ##          SNP effect_allele.exposure other_allele.exposure effect_allele.outcome
    ## 1 rs10471645                      T                     C                     T
    ## 2 rs10472291                      C                     A                     C
    ## 3  rs4871180                      C                     T                     C
    ## 4  rs7624168                      A                     G                     A
    ## 5  rs8074740                      G                     A                     G
    ##   other_allele.outcome beta.exposure beta.outcome eaf.exposure eaf.outcome
    ## 1                    C    0.00487399      -0.0121     0.165686      0.1617
    ## 2                    A   -0.00378453       0.0092     0.666789      0.6632
    ## 3                    T   -0.00382132      -0.0045     0.754276      0.7367
    ## 4                    G   -0.00368970      -0.0023     0.224514      0.2036
    ## 5                    A   -0.00365379      -0.0056     0.677082      0.6640
    ## 
    ## [[2]]
    ##         SNP effect_allele.exposure other_allele.exposure effect_allele.outcome
    ## 1 rs1437405                      T                     C                     T
    ## 2 rs2957297                      C                     A                     C
    ##   other_allele.outcome beta.exposure beta.outcome eaf.exposure eaf.outcome
    ## 1                    C   -0.00105020       0.0144     0.711174      0.6938
    ## 2                    A    0.00126796      -0.0421     0.814293      0.8059

``` r
# remove HP outliers
DIV_HEM_dat_list_adj <- lapply(names(DIV_HEM_dat_list), function(i){
    DIV_HEM_dat_list[[i]] %>% 
    filter(
      ! row_number() %in% DIV_HEM_presso_list[[i]][[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
      ) 
})
names(DIV_HEM_dat_list_adj) <- names(DIV_HEM_dat_list)

# re-test for horizontal pleiotropy
lapply(DIV_HEM_dat_list_adj, mr_pleiotropy_test) %>% 
  bind_rows() %>% select(-1:-2)
```

    ##              outcome                exposure egger_intercept          se
    ## 1 HEM (Zheng et al.) DIV (Schafmayer et al.)    -0.006293818 0.005937565
    ## 2 HEM (Zheng et al.)   DIV (Dönertaş et al.)     0.004352495 0.002714954
    ##        pval
    ## 1 0.2953495
    ## 2 0.1165773

The results show that the selected variants are not likely to be
horizontally pleiotropic (P \< 0.05).

 

#### Perform MR analyses

As previously, using the retained IVs, we perform MR by utilizing
default methods of the `TwoSampleMR::mr()` function.

``` r
# perform MR
DIV_HEM_list_res <- lapply(
  DIV_HEM_dat_list_adj, mr
  )

# odds ratio
DIV_HEM_list_res <- lapply(DIV_HEM_list_res, function(x){
  generate_odds_ratios(x)
})

# results
DIV_HEM_list_res %>% bind_rows() %>% 
  select(-1:-2)
```

    ##               outcome                exposure                    method nsnp
    ## 1  HEM (Zheng et al.) DIV (Schafmayer et al.)                  MR Egger   43
    ## 2  HEM (Zheng et al.) DIV (Schafmayer et al.)           Weighted median   43
    ## 3  HEM (Zheng et al.) DIV (Schafmayer et al.) Inverse variance weighted   43
    ## 4  HEM (Zheng et al.) DIV (Schafmayer et al.)               Simple mode   43
    ## 5  HEM (Zheng et al.) DIV (Schafmayer et al.)             Weighted mode   43
    ## 6  HEM (Zheng et al.)   DIV (Dönertaş et al.)                  MR Egger   43
    ## 7  HEM (Zheng et al.)   DIV (Dönertaş et al.)           Weighted median   43
    ## 8  HEM (Zheng et al.)   DIV (Dönertaş et al.) Inverse variance weighted   43
    ## 9  HEM (Zheng et al.)   DIV (Dönertaş et al.)               Simple mode   43
    ## 10 HEM (Zheng et al.)   DIV (Dönertaş et al.)             Weighted mode   43
    ##             b        se         pval       lo_ci    up_ci         or   or_lci95
    ## 1   2.9042916 1.3276100 3.445506e-02  0.30217607 5.506407 18.2523100 1.35279939
    ## 2   0.7513246 0.2763797 6.558833e-03  0.20962049 1.293029  2.1198062 1.23320996
    ## 3   1.5707530 0.4247044 2.169101e-04  0.73833235 2.403174  4.8102690 2.09244314
    ## 4   0.6704240 0.4925173 1.807032e-01 -0.29490977 1.635758  1.9550662 0.74459877
    ## 5   0.7933646 0.3931688 5.002264e-02  0.02275385 1.563975  2.2108225 1.02301469
    ## 6  -0.3638746 2.0459310 8.597144e-01 -4.37389933 3.646150  0.6949784 0.01260201
    ## 7   2.2847784 1.0331958 2.701022e-02  0.25971458 4.309842  9.8235091 1.29655997
    ## 8   2.7346302 0.6834722 6.305273e-05  1.39502476 4.074236 15.4040467 4.03507448
    ## 9   0.8350969 2.1067664 6.938256e-01 -3.29416528 4.964359  2.3050374 0.03709900
    ## 10  2.2689945 1.6454109 1.751998e-01 -0.95601096 5.494000  9.6696729 0.38442331
    ##      or_uci95
    ## 1  246.264763
    ## 2    3.643806
    ## 3   11.058216
    ## 4    5.133347
    ## 5    4.777777
    ## 6   38.326831
    ## 7   74.428745
    ## 8   58.805521
    ## 9  143.216729
    ## 10 243.228161

``` r
# plot results
DIV_HEM_list_plot <- lapply(seq(DIV_HEM_list_res), function(i){
  mr_scatter_plot(DIV_HEM_list_res[[i]], DIV_HEM_dat_list_adj[[i]])
})
names(DIV_HEM_list_plot) <- names(DIV_HEM_list_res)

lapply(DIV_HEM_list_plot, function(x){
  x[[1]]
})
```

    ## $`DIV (Schafmayer et al.)`

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_refinement_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

    ## 
    ## $`DIV (Dönertaş et al.)`

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_refinement_files/figure-gfm/unnamed-chunk-20-2.png)<!-- -->

In both of the datasets, we observe statistically significant
associations ($P_{IWV}<0.05$; $bxy>0$) of DIV (exposure) on HEM
(outcome).

 

#### Perform sensitivity analyses

Here we perform leave-one-out analysis to test if association was not
driven by a single variant.

``` r
# leave one out analysis
DIV_HEM_list_loo <- lapply(
  DIV_HEM_dat_list_adj, mr_leaveoneout
  )

# plot results
lapply(DIV_HEM_list_loo, function(x){
  p <- mr_leaveoneout_plot(x)
  p[[1]] + coord_flip() +
  theme(
    axis.text.x = element_text(
      angle = 90, vjust=0.5, hjust = 1, size = 6
    ))
  })
```

    ## $`DIV (Schafmayer et al.)`

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_refinement_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

    ## 
    ## $`DIV (Dönertaş et al.)`

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_refinement_files/figure-gfm/unnamed-chunk-21-2.png)<!-- -->

Sensitivity analyses show that there is no single variant/s which are
driving the association.

 

### IBS (exposure) and HEM (outcome)

 

##### Obtain IVs from GWAS of exposure (IBS)

Here we obtain SNPs from the two IBS GWAS summary statistics using
clumping procedure with refined thresholds:

``` r
# get effects of instruments on outcome
IBS_exp_dat_list <- lapply(IBS_file_list, function(x){
  read_exposure_data(
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
    ) %>%
    filter(
      pval.exposure < 1e-05
      )
  })

# add exposure name
for(i in seq(IBS_exp_dat_list)){
  IBS_exp_dat_list[[i]] <- IBS_exp_dat_list[[i]] %>% 
    mutate(
      exposure = names(IBS_exp_dat_list[i])
    )
}

# clump instruments
IBS_exp_dat_list <- lapply(
  IBS_exp_dat_list,
  clump_data,
  clump_r2 = 0.001,
  clump_kb = 1000
  )

# number of SNPs
lapply(IBS_exp_dat_list, nrow)
```

    ## $`IBS (Dönertaş et al.)`
    ## [1] 36
    ## 
    ## $`IBS (Jiang et al.)`
    ## [1] 33

``` r
# get effects of instruments on outcome
HEM_out_dat_list <- lapply(IBS_exp_dat_list, function(x){
  extract_outcome_data(
    snps = x$SNP,
    outcomes = "ebi-a-GCST90014033"
    ) %>% 
    mutate(
      outcome = "HEM (Zheng et al.)"
      )
  })

# check if any SNPs are missing
missing <- lapply(seq(HEM_out_dat_list), function(i){
    IBS_exp_dat_list[[i]] %>% filter(!SNP %in% HEM_out_dat_list[[i]]$SNP) %>% 
    pull(SNP)
  })
names(missing) <- names(HEM_out_dat_list)
missing
```

    ## $`IBS (Dönertaş et al.)`
    ## [1] "rs532907398" "rs532448066" "rs571603499" "rs143829612"
    ## 
    ## $`IBS (Jiang et al.)`
    ## [1] "rs7104964"

 

#### Harmonize exposure and outcome data

As previously, we harmonize SNPs in exposure and outcome datasets.

``` r
# harmonize the exposure and outcome data
IBS_HEM_dat_list <- lapply(seq(IBS_exp_dat_list), function(i){
  harmonise_data(IBS_exp_dat_list[[i]], HEM_out_dat_list[[i]])
})
names(IBS_HEM_dat_list) <- names(HEM_out_dat_list)

# number of SNPs retained for MR analysis
lapply(IBS_HEM_dat_list, function(x){
  x %>% filter(mr_keep == TRUE) %>% 
    nrow()
})
```

    ## $`IBS (Dönertaş et al.)`
    ## [1] 31
    ## 
    ## $`IBS (Jiang et al.)`
    ## [1] 32

 

#### Test for horizontal pleiotropy

As before, we use MR Egger intercept test to evaluate if selected IVs
are not horizontally pleiotropic.

``` r
# test for horizontal pleiotropy
lapply(IBS_HEM_dat_list, mr_pleiotropy_test) %>% 
  bind_rows() %>% select(-1:-2)
```

    ##              outcome              exposure egger_intercept          se
    ## 1 HEM (Zheng et al.) IBS (Dönertaş et al.)    0.0004588635 0.002649512
    ## 2 HEM (Zheng et al.)    IBS (Jiang et al.)   -0.0009797440 0.002535652
    ##        pval
    ## 1 0.8637069
    ## 2 0.7019368

 

Then again, we use MR-PRESSO to detect and remove outliers and then
again run MR Egger intercept test to evaluate horizontal pleiotropy.

``` r
# set random seed
set.seed(0)

# run mr-presso to perform correction of
# horizontal pleiotropy via outlier removal
IBS_HEM_presso_list <- lapply(
  IBS_HEM_dat_list,
  run_mr_presso,
  NbDistribution = 2000,
  SignifThreshold = 0.1
  )

# global test results
lapply(IBS_HEM_presso_list, function(x){
  x[[1]]$`MR-PRESSO results`$`Global Test`
})
```

    ## $`IBS (Dönertaş et al.)`
    ## $`IBS (Dönertaş et al.)`$RSSobs
    ## [1] 31.62585
    ## 
    ## $`IBS (Dönertaş et al.)`$Pvalue
    ## [1] 0.5105
    ## 
    ## 
    ## $`IBS (Jiang et al.)`
    ## $`IBS (Jiang et al.)`$RSSobs
    ## [1] 28.99957
    ## 
    ## $`IBS (Jiang et al.)`$Pvalue
    ## [1] 0.6655

Both methods show that the selected variants are not likely to be
horizontally pleiotropic (P \< 0.05).

 

#### Perform MR analyses

As previously, using the retained IVs, we perform MR by utilizing
default methods of the `TwoSampleMR::mr()` function.

``` r
# perform MR
IBS_HEM_list_res <- lapply(
  IBS_HEM_dat_list, mr
  )

# add odds ratio
IBS_HEM_list_res <- lapply(IBS_HEM_list_res, function(x){
  generate_odds_ratios(x)
})

# results
IBS_HEM_list_res %>% bind_rows() %>% 
  select(-1:-2)
```

    ##               outcome              exposure                    method nsnp
    ## 1  HEM (Zheng et al.) IBS (Dönertaş et al.)                  MR Egger   31
    ## 2  HEM (Zheng et al.) IBS (Dönertaş et al.)           Weighted median   31
    ## 3  HEM (Zheng et al.) IBS (Dönertaş et al.) Inverse variance weighted   31
    ## 4  HEM (Zheng et al.) IBS (Dönertaş et al.)               Simple mode   31
    ## 5  HEM (Zheng et al.) IBS (Dönertaş et al.)             Weighted mode   31
    ## 6  HEM (Zheng et al.)    IBS (Jiang et al.)                  MR Egger   32
    ## 7  HEM (Zheng et al.)    IBS (Jiang et al.)           Weighted median   32
    ## 8  HEM (Zheng et al.)    IBS (Jiang et al.) Inverse variance weighted   32
    ## 9  HEM (Zheng et al.)    IBS (Jiang et al.)               Simple mode   32
    ## 10 HEM (Zheng et al.)    IBS (Jiang et al.)             Weighted mode   32
    ##               b          se      pval        lo_ci      up_ci        or
    ## 1   0.513742854 1.424815969 0.7210357 -2.278896446 3.30638215 1.6715358
    ## 2   0.499408979 0.809688824 0.5373724 -1.087581116 2.08639907 1.6477471
    ## 3   0.739824021 0.565964598 0.1911477 -0.369466592 1.84911463 2.0955667
    ## 4   0.352698402 1.457564862 0.8104439 -2.504128728 3.20952553 1.4229019
    ## 5  -0.067803335 1.339859442 0.9599759 -2.693927842 2.55832117 0.9344442
    ## 6   0.008700773 0.009271347 0.3555036 -0.009471068 0.02687261 1.0087387
    ## 7   0.008099596 0.005903330 0.1700515 -0.003470931 0.01967012 1.0081325
    ## 8   0.005515561 0.004242810 0.1936084 -0.002800347 0.01383147 1.0055308
    ## 9   0.011340231 0.011362408 0.3259831 -0.010930088 0.03361055 1.0114048
    ## 10  0.011340231 0.009626009 0.2477285 -0.007526746 0.03020721 1.0114048
    ##      or_lci95  or_uci95
    ## 1  0.10239715 27.286229
    ## 2  0.33703075  8.055854
    ## 3  0.69110287  6.354191
    ## 4  0.08174679 24.767332
    ## 5  0.06761484 12.914119
    ## 6  0.99057364  1.027237
    ## 7  0.99653509  1.019865
    ## 8  0.99720357  1.013928
    ## 9  0.98912943  1.034182
    ## 10 0.99250151  1.030668

``` r
# plot results
IBS_HEM_list_plot <- lapply(seq(IBS_HEM_list_res), function(i){
  mr_scatter_plot(IBS_HEM_list_res[[i]], IBS_HEM_dat_list[[i]])
})
names(IBS_HEM_list_plot) <- names(IBS_HEM_list_res)

lapply(IBS_HEM_list_plot, function(x){
  x[[1]]
})
```

    ## $`IBS (Dönertaş et al.)`

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_refinement_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

    ## 
    ## $`IBS (Jiang et al.)`

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_refinement_files/figure-gfm/unnamed-chunk-26-2.png)<!-- -->

In both of the datasets, we do not observe any statistically significant
associations ($P_{IWV}<0.05$) of IBS (exposure) on HEM (outcome).
However, in both of the datasets $bxy_{IWV}$ estimate is consistently
positive ($bxy>0$).

 

#### Perform sensitivity analyses

Here we perform leave-one-out analysis to test if association was not
driven by a single variant.

``` r
# leave one out analysis
IBS_HEM_list_loo <- lapply(
  IBS_HEM_dat_list, mr_leaveoneout
  )

# plot results
lapply(IBS_HEM_list_loo, function(x){
  p <- mr_leaveoneout_plot(x)
  p[[1]] + coord_flip() +
  theme(
    axis.text.x = element_text(
      angle = 90, vjust=0.5, hjust = 1, size = 6
    ))
  })
```

    ## $`IBS (Dönertaş et al.)`

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_refinement_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

    ## 
    ## $`IBS (Jiang et al.)`

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_refinement_files/figure-gfm/unnamed-chunk-27-2.png)<!-- -->

 

In conclusion, by refining the parameters of analysis by Zhu *et al.*
study, we could not replicate their findings. Moreover, by performing
reverse MR analysis we show significant causal associations of DIV on
HEM development.

 

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
    ## [16] digest_0.6.29       rvest_1.0.3         MRPRESSO_1.0       
    ## [19] colorspace_2.0-3    R.oo_1.25.0         htmltools_0.5.5    
    ## [22] Matrix_1.5-1        plyr_1.8.7          pkgconfig_2.0.3    
    ## [25] broom_1.0.1         haven_2.5.1         scales_1.2.1       
    ## [28] tzdb_0.3.0          timechange_0.2.0    googledrive_2.0.0  
    ## [31] farver_2.1.1        generics_0.1.3      ellipsis_0.3.2     
    ## [34] withr_2.5.0         cli_3.4.1           survival_3.4-0     
    ## [37] magrittr_2.0.3      crayon_1.5.1        readxl_1.4.1       
    ## [40] ieugwasr_0.1.5      evaluate_0.16       R.methodsS3_1.8.2  
    ## [43] fs_1.5.2            fansi_1.0.3         xml2_1.3.3         
    ## [46] tools_4.2.1         data.table_1.14.2   hms_1.1.2          
    ## [49] gargle_1.2.1        lifecycle_1.0.3     munsell_0.5.0      
    ## [52] reprex_2.0.2        glmnet_4.1-4        mr.raps_0.2        
    ## [55] compiler_4.2.1      rlang_1.1.1         grid_4.2.1         
    ## [58] iterators_1.0.14    rstudioapi_0.14     labeling_0.4.2     
    ## [61] rmarkdown_2.16      gtable_0.3.1        codetools_0.2-18   
    ## [64] DBI_1.1.3           curl_4.3.2          R6_2.5.1           
    ## [67] lubridate_1.9.2     knitr_1.40          fastmap_1.1.0      
    ## [70] utf8_1.2.2          nortest_1.0-4       shape_1.4.6        
    ## [73] stringi_1.7.8       Rcpp_1.0.9          vctrs_0.6.3        
    ## [76] dbplyr_2.2.1        tidyselect_1.2.0    xfun_0.33
