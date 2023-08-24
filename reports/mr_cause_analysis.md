CAUSE MR analysis
================

 

## GWAS summary statistics data

- HEM ([Zheng *et al.*](https://pubmed.ncbi.nlm.nih.gov/33888516/)) -
  GCST90014033
- DIV ([Schafmayer *et
  al.*](https://pubmed.ncbi.nlm.nih.gov/30661054/)) - GCST008105
- IBS ([Eijsbouts *et
  al.*](https://pubmed.ncbi.nlm.nih.gov/34741163/)) - GCST90016564

 

``` r
# libraries
library(cause)
library(ieugwasr)
library(tidyverse)
```

### Download files

``` r
##
## download files
##

ebi_ftp <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"

# DIV - GCST008105
# https://pubmed.ncbi.nlm.nih.gov/30661054/
download.file( 
  url = paste0(ebi_ftp,"GCST008001-GCST009000/GCST008105/GWAS_summary_1-23.dosages.maf_0.01.info_0.4.txt.gz"),
  destfile = "../data/GWAS_summary_1-23.dosages.maf_0.01.info_0.4.txt.gz",
  method = "curl"
)

# IBS - GCST90016564 
# https://pubmed.ncbi.nlm.nih.gov/34741163/
download.file( 
  url = paste0(ebi_ftp,"GCST90016001-GCST90017000/GCST90016564/GCST90016564_buildGRCh37.tsv"),
  destfile = "../data/GCST90016564_buildGRCh37.tsv",
  method = "curl"
)
# compress
R.utils::gzip(
  "../data/GCST90016564_buildGRCh37.tsv",
  overwrite=TRUE
)

# HEM - GCST90014033
# https://pubmed.ncbi.nlm.nih.gov/33888516/
download.file( 
  url = paste0(
    ebi_ftp,"GCST90014001-GCST90015000/GCST90014033/harmonised/",
    "33888516-GCST90014033-EFO_0009552.h.tsv.gz"
    ),
  destfile = "../data/GCST90014033_buildGRCh37_harmonized.tsv.gz",
  method = "curl"
)
```

 

## CAUSE MR analysis

 

### Instrument variable selection

For the instrument variable (IVs) selection (LD clumping) thresholds:

1)  p-value = 1e<sup>-03</sup>;
2)  LD - r<sup>2</sup> = 0.001 (1000 Genomes EUR);

 

### Read data

``` r
# HEM data
HEM <- read_tsv(
  "../data/GCST90014033_buildGRCh37_harmonized.tsv.gz"
)
```

    ## Rows: 8469409 Columns: 24
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (8): hm_variant_id, hm_rsid, hm_chrom, hm_other_allele, hm_effect_allel...
    ## dbl (10): hm_pos, hm_beta, hm_effect_allele_frequency, hm_code, chromosome, ...
    ## lgl  (6): hm_odds_ratio, hm_ci_lower, hm_ci_upper, odds_ratio, ci_lower, ci_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# DIV data
DIV <- read_delim(
  "../data/GWAS_summary_1-23.dosages.maf_0.01.info_0.4.txt.gz",
)
```

    ## New names:
    ## Rows: 10113596 Columns: 16
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: " " chr
    ## (3): SNP, ALLELE1, ALLELE0 dbl (7): CHR, BP, A1FREQ, INFO, BETA, SE, P lgl (6):
    ## ...4, ...9, ...10, ...13, ...14, ...15
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...4`
    ## • `` -> `...9`
    ## • `` -> `...10`
    ## • `` -> `...13`
    ## • `` -> `...14`
    ## • `` -> `...15`

``` r
# IBS data
IBS <- read_tsv(
  "../data/GCST90016564_buildGRCh37.tsv.gz"
)
```

    ## Rows: 9885498 Columns: 11
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): variant_id, effect_allele, other_allele
    ## dbl (8): p_value, chromosome, base_pair_location, effect_allele_frequency, b...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

 

### HEM (exposure) -\> IBS (outcome)

``` r
# harmonize data
HEM_IBS <- gwas_merge(
  HEM, IBS,
  snp_name_cols = c("hm_rsid", "variant_id"),
  beta_hat_cols = c("beta", "beta"), 
  se_cols = c("standard_error", "standard_error"), 
  A1_cols = c("effect_allele", "effect_allele"), 
  A2_cols = c("other_allele", "other_allele"),
  pval_cols = c("p_value", "p_value")
  )
```

    ## Formatting X1
    ## There are  8469409  variants.
    ## Removing  5  duplicated variants leaving  8469399 variants.
    ## Removing  1  variants with illegal alleles leaving  8469399 variants.
    ## Removed  1277502  variants with ambiguous strand.
    ## Flipping strand and effect allele so A1 is always A
    ## Returning  7191897  variants.
    ## Formatting X2
    ## There are  9885498  variants.
    ## Removing  168  duplicated variants leaving  9885162 variants.
    ## Removing  1  variants with illegal alleles leaving  9233338 variants.
    ## Removed  1442205  variants with ambiguous strand.
    ## Flipping strand and effect allele so A1 is always A
    ## Returning  7791133  variants.
    ## After merging and removing variants with inconsistent alleles,  there are  6698725  variants that are present in both studies and can be used with CAUSE.

``` r
## calculate nuisance parameters

# set random seed
set.seed(0)

HEM_IBS_varlist <- with(
  HEM_IBS,
  sample(snp, size=1000000, replace=FALSE)
  )

HEM_IBS_params <- est_cause_params(
  HEM_IBS,
  HEM_IBS_varlist
  )
```

    ## Estimating CAUSE parameters with  1000000  variants.
    ## 1 0.1227584 
    ## 2 0.0002046286 
    ## 3 1.737312e-06 
    ## 4 2.814823e-08

``` r
# LD climping
HEM_IBS_clump <- HEM_IBS %>%
  rename(rsid = "snp", pval = p1) %>%
  ld_clump(
    dat = .,
    clump_r2 = 0.001,
    clump_p = 1e-3,
    plink_bin = genetics.binaRies::get_plink_binary(),
    pop = "EUR"
    )
```

    ## Please look at vignettes for options on running this locally if you need to run many instances of this command.

    ## Clumping 4axZpJ, 6698725 variants, using EUR population reference

    ## Removing 6697981 of 6698725 variants due to LD with other variants or absence from LD reference panel

``` r
# fit CAUSE
HEM_IBS_res <- cause(
  X=HEM_IBS,
  variants = HEM_IBS_clump$rsid,
  param_ests = HEM_IBS_params
  )
```

    ## Estimating CAUSE posteriors using  744  variants.
    ## Fitting confounder only model.
    ## Setting ranges
    ## Refining grid
    ## Fitting causal model.
    ## Setting ranges
    ## Refining grid

``` r
## Pareto k diagnostics

# shared model
HEM_IBS_res$loos[[2]]
```

    ## 
    ## Computed from 1000 by 744 log-likelihood matrix
    ## 
    ##          Estimate    SE
    ## elpd_loo   1056.7  64.6
    ## p_loo         0.9   0.2
    ## looic     -2113.5 129.1
    ## ------
    ## Monte Carlo SE of elpd_loo is 0.0.
    ## 
    ## All Pareto k estimates are good (k < 0.5).
    ## See help('pareto-k-diagnostic') for details.

``` r
# causal model
HEM_IBS_res$loos[[3]]
```

    ## 
    ## Computed from 1000 by 744 log-likelihood matrix
    ## 
    ##          Estimate    SE
    ## elpd_loo   1058.3  64.7
    ## p_loo         1.1   0.2
    ## looic     -2116.5 129.5
    ## ------
    ## Monte Carlo SE of elpd_loo is 0.0.
    ## 
    ## All Pareto k estimates are good (k < 0.5).
    ## See help('pareto-k-diagnostic') for details.

``` r
## results

# expected log pointwise posterior density
HEM_IBS_res$elpd %>% 
  mutate(
    pval = pnorm(z, lower.tail=TRUE)
  )
```

    ##    model1  model2 delta_elpd se_delta_elpd          z      pval
    ## 1    null sharing -0.6784844      1.039004 -0.6530144 0.2568735
    ## 2    null  causal -2.2108448      2.532536 -0.8729768 0.1913379
    ## 3 sharing  causal -1.5323604      1.630972 -0.9395384 0.1737272

``` r
# summary
summary(HEM_IBS_res, ci_size=0.95)
```

    ## p-value testing that causal model is a better fit:  0.17 
    ## Posterior medians and  95 % credible intervals:
    ##      model     gamma            eta                  q               
    ## [1,] "Sharing" NA               "0.33 (-0.24, 1.21)" "0.08 (0, 0.33)"
    ## [2,] "Causal"  "0.07 (0, 0.14)" "0.08 (-0.84, 1.21)" "0.04 (0, 0.25)"

``` r
# visualize
plot(HEM_IBS_res, type="data")
```

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_cause_analysis_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

 

### HEM (exposure) -\> DIV (outcome)

``` r
# harmonize data
HEM_DIV <- gwas_merge(
  HEM, DIV,
  snp_name_cols = c("hm_rsid", "SNP"),
  beta_hat_cols = c("beta", "BETA"), 
  se_cols = c("standard_error", "SE"), 
  A1_cols = c("effect_allele", "ALLELE1"), 
  A2_cols = c("other_allele", "ALLELE0"),
  pval_cols = c("p_value", "P")
  )
```

    ## Formatting X1
    ## There are  8469409  variants.
    ## Removing  5  duplicated variants leaving  8469399 variants.
    ## Removing  1  variants with illegal alleles leaving  8469399 variants.
    ## Removed  1277502  variants with ambiguous strand.
    ## Flipping strand and effect allele so A1 is always A
    ## Returning  7191897  variants.
    ## Formatting X2
    ## There are  10113596  variants.
    ## Removing  12123  duplicated variants leaving  10089669 variants.
    ## Removing  1  variants with illegal alleles leaving  8900692 variants.
    ## Removed  1392521  variants with ambiguous strand.
    ## Flipping strand and effect allele so A1 is always A
    ## Returning  7508171  variants.
    ## After merging and removing variants with inconsistent alleles,  there are  6658231  variants that are present in both studies and can be used with CAUSE.

``` r
## calculate nuisance parameters

# set random seed
set.seed(0)

HEM_DIV_varlist <- with(
  HEM_DIV,
  sample(snp, size=1000000, replace=FALSE)
  )

HEM_DIV_params <- est_cause_params(
  HEM_DIV,
  HEM_DIV_varlist
  )
```

    ## Estimating CAUSE parameters with  1000000  variants.
    ## 1 0.143528 
    ## 2 0.0004749966 
    ## 3 1.147797e-05 
    ## 4 2.778226e-07 
    ## 5 2.210904e-08

``` r
# LD pruning
HEM_DIV_clump <- HEM_DIV %>%
  rename(rsid = "snp", pval = p1) %>%
  ld_clump(
    dat = .,
    clump_r2 = 0.001,
    clump_p = 1e-3,
    plink_bin = genetics.binaRies::get_plink_binary(),
    pop = "EUR"
    )
```

    ## Please look at vignettes for options on running this locally if you need to run many instances of this command.

    ## Clumping ykKqZ2, 6658231 variants, using EUR population reference

    ## Removing 6657491 of 6658231 variants due to LD with other variants or absence from LD reference panel

``` r
# fit CAUSE
HEM_DIV_res <- cause(
  X=HEM_DIV,
  variants = HEM_DIV_clump$rsid,
  param_ests = HEM_DIV_params
  )
```

    ## Estimating CAUSE posteriors using  740  variants.
    ## Fitting confounder only model.
    ## Setting ranges
    ## Refining grid
    ## Fitting causal model.
    ## Setting ranges
    ## Refining grid

``` r
## Pareto k diagnostics

# shared model
HEM_DIV_res$loos[[2]]
```

    ## 
    ## Computed from 1000 by 740 log-likelihood matrix
    ## 
    ##          Estimate    SE
    ## elpd_loo   2911.0  63.7
    ## p_loo         0.6   0.1
    ## looic     -5822.0 127.3
    ## ------
    ## Monte Carlo SE of elpd_loo is 0.0.
    ## 
    ## All Pareto k estimates are good (k < 0.5).
    ## See help('pareto-k-diagnostic') for details.

``` r
# causal model
HEM_DIV_res$loos[[3]]
```

    ## 
    ## Computed from 1000 by 740 log-likelihood matrix
    ## 
    ##          Estimate    SE
    ## elpd_loo   2913.3  63.7
    ## p_loo         0.8   0.1
    ## looic     -5826.6 127.4
    ## ------
    ## Monte Carlo SE of elpd_loo is 0.0.
    ## 
    ## All Pareto k estimates are good (k < 0.5).
    ## See help('pareto-k-diagnostic') for details.

``` r
## results

# expected log pointwise posterior density
HEM_DIV_res$elpd %>% 
  mutate(
    pval = pnorm(z, lower.tail=TRUE)
  )
```

    ##    model1  model2 delta_elpd se_delta_elpd          z       pval
    ## 1    null sharing -0.3832443     0.7463914 -0.5134629 0.30381378
    ## 2    null  causal -2.6559734     2.4486169 -1.0846831 0.13903102
    ## 3 sharing  causal -2.2727291     1.7588640 -1.2921574 0.09815133

``` r
# summary
summary(HEM_DIV_res, ci_size=0.95)
```

    ## p-value testing that causal model is a better fit:  0.098 
    ## Posterior medians and  95 % credible intervals:
    ##      model     gamma            eta                  q               
    ## [1,] "Sharing" NA               "0.03 (-0.03, 0.12)" "0.07 (0, 0.34)"
    ## [2,] "Causal"  "0.01 (0, 0.01)" "0 (-0.08, 0.12)"    "0.03 (0, 0.24)"

``` r
# visualize
plot(HEM_DIV_res, type="data")
```

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_cause_analysis_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

 

## Reverse MR analysis

### IBS (exposure) -\> HEM (outcome)

``` r
# harmonize data
IBS_HEM <- gwas_merge(
  IBS, HEM,
  snp_name_cols = c("variant_id", "hm_rsid"),
  beta_hat_cols = c("beta", "beta"), 
  se_cols = c("standard_error", "standard_error"), 
  A1_cols = c("effect_allele", "effect_allele"), 
  A2_cols = c("other_allele", "other_allele"),
  pval_cols = c("p_value", "p_value")
  )
```

    ## Formatting X1
    ## There are  9885498  variants.
    ## Removing  168  duplicated variants leaving  9885162 variants.
    ## Removing  1  variants with illegal alleles leaving  9233338 variants.
    ## Removed  1442205  variants with ambiguous strand.
    ## Flipping strand and effect allele so A1 is always A
    ## Returning  7791133  variants.
    ## Formatting X2
    ## There are  8469409  variants.
    ## Removing  5  duplicated variants leaving  8469399 variants.
    ## Removing  1  variants with illegal alleles leaving  8469399 variants.
    ## Removed  1277502  variants with ambiguous strand.
    ## Flipping strand and effect allele so A1 is always A
    ## Returning  7191897  variants.
    ## After merging and removing variants with inconsistent alleles,  there are  6698725  variants that are present in both studies and can be used with CAUSE.

``` r
## calculate nuisance parameters

# set random seed
set.seed(0)

IBS_HEM_varlist <- with(
  IBS_HEM,
  sample(snp, size=1000000, replace=FALSE)
  )

IBS_HEM_params <- est_cause_params(
  IBS_HEM,
  IBS_HEM_varlist
  )
```

    ## Estimating CAUSE parameters with  1000000  variants.
    ## 1 0.1072283 
    ## 2 0.01332328 
    ## 3 7.310684e-06 
    ## 4 7.983068e-08

``` r
# LD pruning
IBS_HEM_clump <- IBS_HEM %>%
  rename(rsid = "snp", pval = p1) %>%
  ld_clump(
    dat = .,
    clump_r2 = 0.001,
    clump_p = 1e-3,
    plink_bin = genetics.binaRies::get_plink_binary(),
    pop = "EUR"
    )
```

    ## Please look at vignettes for options on running this locally if you need to run many instances of this command.

    ## Clumping 4axZpJ, 6698725 variants, using EUR population reference

    ## Removing 6698117 of 6698725 variants due to LD with other variants or absence from LD reference panel

``` r
# fit CAUSE
IBS_HEM_res <- cause(
  X=IBS_HEM,
  variants = IBS_HEM_clump$rsid,
  param_ests = IBS_HEM_params
  )
```

    ## Estimating CAUSE posteriors using  608  variants.
    ## Fitting confounder only model.
    ## Setting ranges
    ## Refining grid
    ## Fitting causal model.
    ## Setting ranges
    ## Refining grid

``` r
## Pareto k diagnostics

# shared model
IBS_HEM_res$loos[[2]]
```

    ## 
    ## Computed from 1000 by 608 log-likelihood matrix
    ## 
    ##          Estimate    SE
    ## elpd_loo     22.2  55.1
    ## p_loo         0.5   0.1
    ## looic       -44.4 110.2
    ## ------
    ## Monte Carlo SE of elpd_loo is 0.0.
    ## 
    ## All Pareto k estimates are good (k < 0.5).
    ## See help('pareto-k-diagnostic') for details.

``` r
# causal model
IBS_HEM_res$loos[[3]]
```

    ## 
    ## Computed from 1000 by 608 log-likelihood matrix
    ## 
    ##          Estimate    SE
    ## elpd_loo     25.7  55.2
    ## p_loo         0.7   0.1
    ## looic       -51.5 110.4
    ## ------
    ## Monte Carlo SE of elpd_loo is 0.0.
    ## 
    ## All Pareto k estimates are good (k < 0.5).
    ## See help('pareto-k-diagnostic') for details.

``` r
## results

# expected log pointwise posterior density
IBS_HEM_res$elpd %>% 
  mutate(
    pval = pnorm(z, lower.tail=TRUE)
  )
```

    ##    model1  model2 delta_elpd se_delta_elpd         z       pval
    ## 1    null sharing -0.9959641     0.7074238 -1.407875 0.07958410
    ## 2    null  causal -4.5264739     2.7882096 -1.623434 0.05224835
    ## 3 sharing  causal -3.5305098     2.0979216 -1.682861 0.04620103

``` r
# summary
summary(IBS_HEM_res, ci_size=0.95)
```

    ## p-value testing that causal model is a better fit:  0.046 
    ## Posterior medians and  95 % credible intervals:
    ##      model     gamma               eta                   q               
    ## [1,] "Sharing" NA                  "0.2 (-0.24, 0.55)"   "0.12 (0, 0.42)"
    ## [2,] "Causal"  "0.11 (0.04, 0.18)" "-0.01 (-0.54, 0.49)" "0.04 (0, 0.26)"

``` r
# visualize
plot(IBS_HEM_res, type="data")
```

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_cause_analysis_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
 

### DIV (exposure) -\> HEM (outcome)

``` r
# harmonize data
DIV_HEM <- gwas_merge(
  DIV, HEM,
  snp_name_cols = rev(c("hm_rsid", "SNP")),
  beta_hat_cols = rev(c("beta", "BETA")), 
  se_cols = rev(c("standard_error", "SE")), 
  A1_cols = rev(c("effect_allele", "ALLELE1")), 
  A2_cols = rev(c("other_allele", "ALLELE0")),
  pval_cols = rev(c("p_value", "P"))
  )
```

    ## Formatting X1
    ## There are  10113596  variants.
    ## Removing  12123  duplicated variants leaving  10089669 variants.
    ## Removing  1  variants with illegal alleles leaving  8900692 variants.
    ## Removed  1392521  variants with ambiguous strand.
    ## Flipping strand and effect allele so A1 is always A
    ## Returning  7508171  variants.
    ## Formatting X2
    ## There are  8469409  variants.
    ## Removing  5  duplicated variants leaving  8469399 variants.
    ## Removing  1  variants with illegal alleles leaving  8469399 variants.
    ## Removed  1277502  variants with ambiguous strand.
    ## Flipping strand and effect allele so A1 is always A
    ## Returning  7191897  variants.
    ## After merging and removing variants with inconsistent alleles,  there are  6658231  variants that are present in both studies and can be used with CAUSE.

``` r
## calculate nuisance parameters

# set random seed
set.seed(0)

DIV_HEM_varlist <- with(
  DIV_HEM,
  sample(snp, size=1000000, replace=FALSE)
  )

DIV_HEM_params <- est_cause_params(
  DIV_HEM,
  DIV_HEM_varlist
  )
```

    ## Estimating CAUSE parameters with  1000000  variants.
    ## 1 0.1448956 
    ## 2 0.0004638866 
    ## 3 1.131384e-05 
    ## 4 2.770741e-07 
    ## 5 1.618336e-08

``` r
# LD pruning
DIV_HEM_clump <- DIV_HEM %>%
  rename(rsid = "snp", pval = p1) %>%
  ld_clump(
    dat = .,
    clump_r2 = 0.001,
    clump_p = 1e-3,
    plink_bin = genetics.binaRies::get_plink_binary(),
    pop = "EUR"
    )
```

    ## Please look at vignettes for options on running this locally if you need to run many instances of this command.

    ## Clumping ykKqZ2, 6658231 variants, using EUR population reference

    ## Removing 6657579 of 6658231 variants due to LD with other variants or absence from LD reference panel

``` r
# fit CAUSE
DIV_HEM_res <- cause(
  X=DIV_HEM,
  variants = DIV_HEM_clump$rsid,
  param_ests = DIV_HEM_params
  )
```

    ## Estimating CAUSE posteriors using  652  variants.
    ## Fitting confounder only model.
    ## Setting ranges
    ## Refining grid
    ## Fitting causal model.
    ## Setting ranges
    ## Refining grid

``` r
## Pareto k diagnostics

# shared model
DIV_HEM_res$loos[[2]]
```

    ## 
    ## Computed from 1000 by 652 log-likelihood matrix
    ## 
    ##          Estimate    SE
    ## elpd_loo   2113.6  58.6
    ## p_loo         0.6   0.1
    ## looic     -4227.3 117.1
    ## ------
    ## Monte Carlo SE of elpd_loo is 0.0.
    ## 
    ## Pareto k diagnostic values:
    ##                          Count Pct.    Min. n_eff
    ## (-Inf, 0.5]   (good)     650   99.7%   974       
    ##  (0.5, 0.7]   (ok)         2    0.3%   1000      
    ##    (0.7, 1]   (bad)        0    0.0%   <NA>      
    ##    (1, Inf)   (very bad)   0    0.0%   <NA>      
    ## 
    ## All Pareto k estimates are ok (k < 0.7).
    ## See help('pareto-k-diagnostic') for details.

``` r
# causal model
DIV_HEM_res$loos[[3]]
```

    ## 
    ## Computed from 1000 by 652 log-likelihood matrix
    ## 
    ##          Estimate    SE
    ## elpd_loo   2118.7  58.3
    ## p_loo         0.7   0.1
    ## looic     -4237.5 116.6
    ## ------
    ## Monte Carlo SE of elpd_loo is 0.0.
    ## 
    ## All Pareto k estimates are good (k < 0.5).
    ## See help('pareto-k-diagnostic') for details.

``` r
## results

# expected log pointwise posterior density
DIV_HEM_res$elpd %>% 
  mutate(
    pval = pnorm(z, lower.tail=TRUE)
  )
```

    ##    model1  model2 delta_elpd se_delta_elpd         z        pval
    ## 1    null sharing  -1.591004     0.9583344 -1.660176 0.048439503
    ## 2    null  causal  -6.687409     3.0556590 -2.188532 0.014315419
    ## 3 sharing  causal  -5.096405     2.1770198 -2.341001 0.009616065

``` r
# summary
summary(DIV_HEM_res, ci_size=0.95)
```

    ## p-value testing that causal model is a better fit:  0.0096 
    ## Posterior medians and  95 % credible intervals:
    ##      model     gamma               eta                  q               
    ## [1,] "Sharing" NA                  "1.89 (-1.31, 6.77)" "0.14 (0, 0.46)"
    ## [2,] "Causal"  "0.94 (0.44, 1.42)" "0.17 (-6.07, 7.62)" "0.03 (0, 0.24)"

``` r
# visualize
plot(DIV_HEM_res, type="data")
```

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_cause_analysis_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->
 

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
    ##  [1] forcats_0.5.2    stringr_1.4.1    dplyr_1.1.2      purrr_0.3.4     
    ##  [5] readr_2.1.2      tidyr_1.2.1      tibble_3.2.1     ggplot2_3.4.2   
    ##  [9] tidyverse_1.3.2  ieugwasr_0.1.5   cause_1.2.0.0335
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] httr_1.4.4          bit64_4.0.5         vroom_1.5.7        
    ##  [4] jsonlite_1.8.0      modelr_0.1.9        RcppParallel_5.1.7 
    ##  [7] assertthat_0.2.1    highr_0.9           mixsqp_0.3-48      
    ## [10] googlesheets4_1.0.1 cellranger_1.1.0    yaml_2.3.5         
    ## [13] pillar_1.9.0        backports_1.4.1     lattice_0.20-45    
    ## [16] glue_1.6.2          digest_0.6.29       rvest_1.0.3        
    ## [19] colorspace_2.0-3    htmltools_0.5.5     Matrix_1.5-1       
    ## [22] pkgconfig_2.0.3     invgamma_1.1        broom_1.0.1        
    ## [25] haven_2.5.1         scales_1.2.1        intervals_0.15.4   
    ## [28] tzdb_0.3.0          timechange_0.2.0    googledrive_2.0.0  
    ## [31] farver_2.1.1        generics_0.1.3      ellipsis_0.3.2     
    ## [34] withr_2.5.0         ashr_2.2-54         cli_3.4.1          
    ## [37] crayon_1.5.1        magrittr_2.0.3      readxl_1.4.1       
    ## [40] evaluate_0.16       fs_1.5.2            fansi_1.0.3        
    ## [43] xml2_1.3.3          truncnorm_1.0-9     tools_4.2.1        
    ## [46] loo_2.6.0           hms_1.1.2           gargle_1.2.1       
    ## [49] lifecycle_1.0.3     matrixStats_0.62.0  munsell_0.5.0      
    ## [52] reprex_2.0.2        irlba_2.3.5.1       compiler_4.2.1     
    ## [55] rlang_1.1.1         grid_4.2.1          rstudioapi_0.14    
    ## [58] labeling_0.4.2      rmarkdown_2.16      gtable_0.3.1       
    ## [61] curl_4.3.2          DBI_1.1.3           R6_2.5.1           
    ## [64] gridExtra_2.3       lubridate_1.9.2     knitr_1.40         
    ## [67] bit_4.0.4           fastmap_1.1.0       utf8_1.2.2         
    ## [70] stringi_1.7.8       parallel_4.2.1      SQUAREM_2021.1     
    ## [73] Rcpp_1.0.9          vctrs_0.6.3         dbplyr_2.2.1       
    ## [76] tidyselect_1.2.0    xfun_0.33
