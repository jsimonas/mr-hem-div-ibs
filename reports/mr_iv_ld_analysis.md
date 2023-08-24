Linkage disequilirium analysis
================

## Introduction

In the study by [Zhu *et
al.*](https://gut.bmj.com/content/early/2023/01/23/gutjnl-2022-329307-0),
the number of IVs used in Mendelian Randomization (MR) analyses were
always exceeding the number of independent loci, which were identified
in the original HEM GWAS study of [Zheng *et
al.*](https://gut.bmj.com/content/70/8/1538) This observation raised a
concern, because using correlating genetic variants (that is, in linkage
disequilibrium) is breaching the core assumptions of the IWV method,
which requires to use independent IVs [Burgess *et
al.*](https://link.springer.com/article/10.1007/s10654-017-0255-x)

Therefore, we have tested if correlation threshold (r<sup>2</sup> \<
0.001) within 10kb window is enough to obtain uncorrelated variants from
HEM GWAS summary statistics data.

 

## Methods

### SNP clumping for LD

Thresholds used by Zhu *et al.*:

1)  Association with exposure (HEM GWAS) - P \< 5e<sup>-10</sup>;
2)  Pruning for linkage disequilibrium (LD) - r<sup>2</sup> \< 0.001
    within 10kb (1000G EUR);
3)  Association with outcome - P \< 0.05.

Commonly used thresholds for LD window:

1)  Clumping for linkage disequilibrium (LD) - r<sup>2</sup> \< 0.001
    within 1 Mb (1000G EUR);
2)  Clumping for linkage disequilibrium (LD) - r<sup>2</sup> \< 0.001
    within 10 Mb (1000G EUR);

 

### Packages

Analyses were performed using R packages
[`LDlinkR`](https://CRAN.R-project.org/package=LDlinkR),
[`TwoSampleMR`](https://mrcieu.github.io/TwoSampleMR/) and
[`tidyverse`](https://www.tidyverse.org).

``` r
library(tidyverse)
library(patchwork)
library(TwoSampleMR)
library(LDlinkR)
```

 

## Analysis and Results

 

##### SNP clumping for LD in the HEM GWAS data

Here we obtain SNPs from HEM GWAS using clumping thresholds as defined
by Zhu *et al.* (10 kb), as well as [commonly used thresholds](ref) (1
mb and 10 mb). To be more precise, only the clumping distance cutoff was
changed in older to evaluate its effect:

``` r
# set distance thresholds (in kilobases)
kb <- c(10, 1e3)

# data from OpenGWAS db
hem_iv_clumping <- lapply(kb, function(x){
  if(x == 10){
    pval = 5e-10
  } else{
    pval = 5e-8
  }
  extract_instruments(
    "ebi-a-GCST90014033",
    p1 = pval,
    r2 = 0.001,
    kb = x
  )
})
names(hem_iv_clumping) <- c(
  "10 Kb and P < 5e-10 cutoff ",
  "1 Mb and P < 5e-8 cutoff"
)

# number of SNPs after clumping
lapply(hem_iv_clumping, nrow)
```

    ## $`10 Kb and P < 5e-10 cutoff `
    ## [1] 318
    ## 
    ## $`1 Mb and P < 5e-8 cutoff`
    ## [1] 101

Here we import DIV and IBS (outcome) GWAS datasets to remove SNPs, which
are associated (P \< 0.05) with these traits.

``` r
# list files to studies
out_file_list <- list(
  "DIV (Schafmayer et al.)" = "../data/GWAS_summary_1-23.dosages.maf_0.01.info_0.4.txt.gz",
  "DIV (Dönertaş et al.)" = "../data/GCST90038682_buildGRCh37.tsv.gz",
  "IBS (Dönertaş et al.)" = "../data/GCST90038626_buildGRCh37.tsv.gz",
  "IBS (Jiang et al.)" = "../data/GCST90044173_buildGRCh37.tsv.gz"
)

# read outcome data
out_dat_list <- lapply(hem_iv_clumping, function(x){
  lapply(seq(out_file_list), function(i){
    if(i==1){
      read_outcome_data(
        snps = x$SNP,
        filename = out_file_list[[i]],
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
        outcome = names(out_file_list)[i]
        )
      } else {
        read_outcome_data(
          snps = x$SNP,
          filename = out_file_list[[i]],
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
        outcome = names(out_file_list)[i]
        )
      }
  })
})

# harmonize data and exclude outcome-associated SNPs (P<0.05)
# as described by Zhu et al.
har_dat_list <- lapply(seq(out_dat_list), function(i){
  lapply(seq(out_dat_list[[i]]), function(j){
    harmonise_data(hem_iv_clumping[[i]], out_dat_list[[i]][[j]]) %>%
      filter(pval.outcome > 0.05 & mr_keep == TRUE)
    }) %>% 
    bind_rows() %>%
    select(SNP, chr.exposure) %>% 
    distinct()
  })
names(har_dat_list) <- names(out_dat_list)

# number of SNPs used as IVs
lapply(har_dat_list, nrow)
```

    ## $`10 Kb and P < 5e-10 cutoff `
    ## [1] 301
    ## 
    ## $`1 Mb and P < 5e-8 cutoff`
    ## [1] 96

``` r
# please provide your LDlink token
# copy the token in the following file:
token <- readLines("../data/LDlink_token.txt")

# split SNPs by chrom
snp_list <- lapply(har_dat_list, function(x){
  y <- split(x, x$chr.exposure)
  # filter chr having only one variant
  idx <- lapply(y, function(i) nrow(i)!=1)
  y <- y[unlist(idx)]
  return(y)
}) 

# retrieve LD
ld_mat_list <- lapply(snp_list, function(i){
  lapply(i, function(x){
    m <- LDmatrix(
      x$SNP,
      pop = "CEU",
      r2d = "r2",
      token = token
      )
    m <- as.matrix(m[,-1])
    m[lower.tri(m, diag = TRUE)] <- NA
    rownames(m) <- colnames(m)
    return(m)
  })
})

# make tabular
ld_mat_tb <- lapply(ld_mat_list, function(i){
  lapply(i, function(x){
    x %>% reshape2::melt(
      varnames = c("SNP1", "SNP2"),
      value.name = "r2"
    )
    }) %>%
  bind_rows(.id = "chr")
  }) %>%
  bind_rows(.id = "dist_cutoff") %>% 
  mutate(
    dist_cutoff = factor(
      dist_cutoff, levels = names(hem_iv_clumping)
    ),
    chr = factor(chr, levels = 1:22)
  )
```

``` r
# make plot
p1 <- ld_mat_tb %>%
  filter(grepl("Kb", dist_cutoff)) %>% 
  ggplot(aes(SNP1, SNP2)) +
  geom_tile(aes(fill = r2)) +
  scale_fill_viridis_c(
    na.value = NA, direction = -1, option = 3, limits = c(0, 1)
    ) +
  theme_bw() +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank()
  ) +
  labs(
    title = expression(10~Kb~distance~cutoff~and~P<5%*%10^-10)
  )


p2 <- ld_mat_tb %>%
  filter(grepl("Mb", dist_cutoff)) %>% 
  mutate(
    SNP1 = factor(SNP1, levels = unique(SNP1)),
    SNP2 = factor(SNP2, levels = unique(SNP2)),
  ) %>% 
  ggplot(aes(SNP1, SNP2)) +
  geom_tile(aes(fill = r2)) +
  scale_fill_viridis_c(
    na.value = NA, direction = -1, option = 3, limits = c(0, 1)
    ) +
  theme_bw() +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank()
  ) +
  labs(
    title = expression(1~Mb~distance~cutoff~and~P<5%*%10^-8)
  )

p1 + p2 + plot_layout(guides='collect') &
  theme(legend.position='bottom')
```

![](/Users/simonasj/Documents/Projects/git/mr-hem-div-ibs/reports/mr_iv_ld_analysis_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Here we list SNP pairs, for which r<sup>2</sup> values are equal to one!

``` r
ld_mat_tb %>% 
  filter(r2 == 1) %>%
  as.data.frame()
```

    ##                     dist_cutoff chr        SNP1        SNP2 r2
    ## 1   10 Kb and P < 5e-10 cutoff    1  rs11578225  rs11581918  1
    ## 2   10 Kb and P < 5e-10 cutoff    1 rs145163454 rs144737447  1
    ## 3   10 Kb and P < 5e-10 cutoff    1 rs145163454   rs2227246  1
    ## 4   10 Kb and P < 5e-10 cutoff    1 rs144737447   rs2227246  1
    ## 5   10 Kb and P < 5e-10 cutoff    1  rs77979353  rs10800428  1
    ## 6   10 Kb and P < 5e-10 cutoff    1 rs145163454   rs1209731  1
    ## 7   10 Kb and P < 5e-10 cutoff    1 rs144737447   rs1209731  1
    ## 8   10 Kb and P < 5e-10 cutoff    1   rs2227246   rs1209731  1
    ## 9   10 Kb and P < 5e-10 cutoff    1   rs1894692      rs6025  1
    ## 10  10 Kb and P < 5e-10 cutoff   12  rs74097857  rs11176001  1
    ## 11  10 Kb and P < 5e-10 cutoff   13   rs1262776   rs1626445  1
    ## 12  10 Kb and P < 5e-10 cutoff   15  rs17293632  rs56062135  1
    ## 13  10 Kb and P < 5e-10 cutoff    2 rs113645544   rs7575552  1
    ## 14  10 Kb and P < 5e-10 cutoff    2  rs57640566  rs79433506  1
    ## 15  10 Kb and P < 5e-10 cutoff    2  rs57640566   rs1113419  1
    ## 16  10 Kb and P < 5e-10 cutoff    2  rs79433506   rs1113419  1
    ## 17  10 Kb and P < 5e-10 cutoff    2  rs57640566 rs113494688  1
    ## 18  10 Kb and P < 5e-10 cutoff    2  rs79433506 rs113494688  1
    ## 19  10 Kb and P < 5e-10 cutoff    2   rs1113419 rs113494688  1
    ## 20  10 Kb and P < 5e-10 cutoff    2  rs57640566 rs113629868  1
    ## 21  10 Kb and P < 5e-10 cutoff    2  rs79433506 rs113629868  1
    ## 22  10 Kb and P < 5e-10 cutoff    2   rs1113419 rs113629868  1
    ## 23  10 Kb and P < 5e-10 cutoff    2 rs113494688 rs113629868  1
    ## 24  10 Kb and P < 5e-10 cutoff    2  rs57640566  rs79464516  1
    ## 25  10 Kb and P < 5e-10 cutoff    2  rs79433506  rs79464516  1
    ## 26  10 Kb and P < 5e-10 cutoff    2   rs1113419  rs79464516  1
    ## 27  10 Kb and P < 5e-10 cutoff    2 rs113494688  rs79464516  1
    ## 28  10 Kb and P < 5e-10 cutoff    2 rs113629868  rs79464516  1
    ## 29  10 Kb and P < 5e-10 cutoff    2   rs4292050   rs6731260  1
    ## 30  10 Kb and P < 5e-10 cutoff    2   rs2124969   rs7567781  1
    ## 31  10 Kb and P < 5e-10 cutoff    2  rs59019419  rs59541669  1
    ## 32  10 Kb and P < 5e-10 cutoff    3   rs2336162   rs2581797  1
    ## 33  10 Kb and P < 5e-10 cutoff    3   rs2336162   rs4519686  1
    ## 34  10 Kb and P < 5e-10 cutoff    3   rs2581797   rs4519686  1
    ## 35  10 Kb and P < 5e-10 cutoff    3   rs2336162   rs9846976  1
    ## 36  10 Kb and P < 5e-10 cutoff    3   rs2581797   rs9846976  1
    ## 37  10 Kb and P < 5e-10 cutoff    3   rs4519686   rs9846976  1
    ## 38  10 Kb and P < 5e-10 cutoff    3   rs2564938   rs1529544  1
    ## 39  10 Kb and P < 5e-10 cutoff    3   rs9847710   rs4687697  1
    ## 40  10 Kb and P < 5e-10 cutoff    3   rs6438003   rs6792493  1
    ## 41  10 Kb and P < 5e-10 cutoff    3  rs56394279   rs3851366  1
    ## 42  10 Kb and P < 5e-10 cutoff    4   rs6533183   rs6839705  1
    ## 43  10 Kb and P < 5e-10 cutoff    4  rs17824374   rs7674065  1
    ## 44  10 Kb and P < 5e-10 cutoff    4   rs7661046   rs4340757  1
    ## 45  10 Kb and P < 5e-10 cutoff    4   rs7661046   rs7378179  1
    ## 46  10 Kb and P < 5e-10 cutoff    4   rs4340757   rs7378179  1
    ## 47  10 Kb and P < 5e-10 cutoff    4   rs7661046  rs10029738  1
    ## 48  10 Kb and P < 5e-10 cutoff    4   rs4340757  rs10029738  1
    ## 49  10 Kb and P < 5e-10 cutoff    4   rs7378179  rs10029738  1
    ## 50  10 Kb and P < 5e-10 cutoff    4   rs7661046  rs12645910  1
    ## 51  10 Kb and P < 5e-10 cutoff    4   rs4340757  rs12645910  1
    ## 52  10 Kb and P < 5e-10 cutoff    4   rs7378179  rs12645910  1
    ## 53  10 Kb and P < 5e-10 cutoff    4  rs10029738  rs12645910  1
    ## 54  10 Kb and P < 5e-10 cutoff    4   rs7661046   rs4376087  1
    ## 55  10 Kb and P < 5e-10 cutoff    4   rs4340757   rs4376087  1
    ## 56  10 Kb and P < 5e-10 cutoff    4   rs7378179   rs4376087  1
    ## 57  10 Kb and P < 5e-10 cutoff    4  rs10029738   rs4376087  1
    ## 58  10 Kb and P < 5e-10 cutoff    4  rs12645910   rs4376087  1
    ## 59  10 Kb and P < 5e-10 cutoff    4   rs7661046  rs17019341  1
    ## 60  10 Kb and P < 5e-10 cutoff    4   rs4340757  rs17019341  1
    ## 61  10 Kb and P < 5e-10 cutoff    4   rs7378179  rs17019341  1
    ## 62  10 Kb and P < 5e-10 cutoff    4  rs10029738  rs17019341  1
    ## 63  10 Kb and P < 5e-10 cutoff    4  rs12645910  rs17019341  1
    ## 64  10 Kb and P < 5e-10 cutoff    4   rs4376087  rs17019341  1
    ## 65  10 Kb and P < 5e-10 cutoff    4   rs7661046  rs72731556  1
    ## 66  10 Kb and P < 5e-10 cutoff    4   rs4340757  rs72731556  1
    ## 67  10 Kb and P < 5e-10 cutoff    4   rs7378179  rs72731556  1
    ## 68  10 Kb and P < 5e-10 cutoff    4  rs10029738  rs72731556  1
    ## 69  10 Kb and P < 5e-10 cutoff    4  rs12645910  rs72731556  1
    ## 70  10 Kb and P < 5e-10 cutoff    4   rs4376087  rs72731556  1
    ## 71  10 Kb and P < 5e-10 cutoff    4  rs17019341  rs72731556  1
    ## 72  10 Kb and P < 5e-10 cutoff    4   rs7661046 rs144044336  1
    ## 73  10 Kb and P < 5e-10 cutoff    4   rs4340757 rs144044336  1
    ## 74  10 Kb and P < 5e-10 cutoff    4   rs7378179 rs144044336  1
    ## 75  10 Kb and P < 5e-10 cutoff    4  rs10029738 rs144044336  1
    ## 76  10 Kb and P < 5e-10 cutoff    4  rs12645910 rs144044336  1
    ## 77  10 Kb and P < 5e-10 cutoff    4   rs4376087 rs144044336  1
    ## 78  10 Kb and P < 5e-10 cutoff    4  rs17019341 rs144044336  1
    ## 79  10 Kb and P < 5e-10 cutoff    4  rs72731556 rs144044336  1
    ## 80  10 Kb and P < 5e-10 cutoff    4  rs62343635   rs2200943  1
    ## 81  10 Kb and P < 5e-10 cutoff    4   rs6845536   rs7663578  1
    ## 82  10 Kb and P < 5e-10 cutoff    5   rs4527146   rs4245972  1
    ## 83  10 Kb and P < 5e-10 cutoff    5   rs4527146  rs72703070  1
    ## 84  10 Kb and P < 5e-10 cutoff    5   rs4245972  rs72703070  1
    ## 85  10 Kb and P < 5e-10 cutoff    5   rs3749618   rs4957080  1
    ## 86  10 Kb and P < 5e-10 cutoff    5   rs3749618 rs113896354  1
    ## 87  10 Kb and P < 5e-10 cutoff    5   rs4957080 rs113896354  1
    ## 88  10 Kb and P < 5e-10 cutoff    5   rs7701852   rs3811910  1
    ## 89  10 Kb and P < 5e-10 cutoff    5   rs6449607  rs62366958  1
    ## 90  10 Kb and P < 5e-10 cutoff    5   rs7701852  rs62368263  1
    ## 91  10 Kb and P < 5e-10 cutoff    5   rs3811910  rs62368263  1
    ## 92  10 Kb and P < 5e-10 cutoff    5   rs6449607  rs62368271  1
    ## 93  10 Kb and P < 5e-10 cutoff    5  rs62366958  rs62368271  1
    ## 94  10 Kb and P < 5e-10 cutoff    5   rs7701852  rs17824230  1
    ## 95  10 Kb and P < 5e-10 cutoff    5   rs3811910  rs17824230  1
    ## 96  10 Kb and P < 5e-10 cutoff    5  rs62368263  rs17824230  1
    ## 97  10 Kb and P < 5e-10 cutoff    6   rs1853148   rs9372480  1
    ## 98  10 Kb and P < 5e-10 cutoff    7   rs4143205  rs10281352  1
    ## 99  10 Kb and P < 5e-10 cutoff    7   rs4143205   rs6462976  1
    ## 100 10 Kb and P < 5e-10 cutoff    7  rs10281352   rs6462976  1
    ## 101 10 Kb and P < 5e-10 cutoff    7   rs4143205   rs6965764  1
    ## 102 10 Kb and P < 5e-10 cutoff    7  rs10281352   rs6965764  1
    ## 103 10 Kb and P < 5e-10 cutoff    7   rs6462976   rs6965764  1
    ## 104 10 Kb and P < 5e-10 cutoff    7  rs17171704   rs3890819  1
    ## 105 10 Kb and P < 5e-10 cutoff    7   rs3757582   rs3801458  1
    ## 106 10 Kb and P < 5e-10 cutoff    7   rs2411046   rs7778418  1
    ## 107 10 Kb and P < 5e-10 cutoff    7  rs10257317  rs10279449  1
    ## 108 10 Kb and P < 5e-10 cutoff    7  rs10257317  rs10271184  1
    ## 109 10 Kb and P < 5e-10 cutoff    7  rs10279449  rs10271184  1
    ## 110 10 Kb and P < 5e-10 cutoff    7   rs7794668  rs10241865  1
    ## 111 10 Kb and P < 5e-10 cutoff    7  rs28856331    rs847648  1
    ## 112 10 Kb and P < 5e-10 cutoff    7  rs28856331   rs9718453  1
    ## 113 10 Kb and P < 5e-10 cutoff    7    rs847648   rs9718453  1
    ## 114 10 Kb and P < 5e-10 cutoff    7  rs28856331  rs11514917  1
    ## 115 10 Kb and P < 5e-10 cutoff    7    rs847648  rs11514917  1
    ## 116 10 Kb and P < 5e-10 cutoff    7   rs9718453  rs11514917  1
    ## 117 10 Kb and P < 5e-10 cutoff    7    rs869332   rs4729866  1
    ## 118 10 Kb and P < 5e-10 cutoff    7  rs10808120  rs12666317  1
    ## 119 10 Kb and P < 5e-10 cutoff    7   rs4729873  rs13221979  1
    ## 120 10 Kb and P < 5e-10 cutoff    7    rs712707    rs712715  1
    ## 121 10 Kb and P < 5e-10 cutoff    7    rs712707    rs712718  1
    ## 122 10 Kb and P < 5e-10 cutoff    7    rs712715    rs712718  1
    ## 123 10 Kb and P < 5e-10 cutoff    7    rs712707    rs712719  1
    ## 124 10 Kb and P < 5e-10 cutoff    7    rs712715    rs712719  1
    ## 125 10 Kb and P < 5e-10 cutoff    7    rs712718    rs712719  1
    ## 126 10 Kb and P < 5e-10 cutoff    7    rs712707  rs12540484  1
    ## 127 10 Kb and P < 5e-10 cutoff    7    rs712715  rs12540484  1
    ## 128 10 Kb and P < 5e-10 cutoff    7    rs712718  rs12540484  1
    ## 129 10 Kb and P < 5e-10 cutoff    7    rs712719  rs12540484  1
    ## 130 10 Kb and P < 5e-10 cutoff    7   rs4731377  rs11772303  1
    ## 131 10 Kb and P < 5e-10 cutoff    7   rs4731377  rs62481385  1
    ## 132 10 Kb and P < 5e-10 cutoff    7  rs11772303  rs62481385  1
    ## 133 10 Kb and P < 5e-10 cutoff    7   rs4731377  rs67338333  1
    ## 134 10 Kb and P < 5e-10 cutoff    7  rs11772303  rs67338333  1
    ## 135 10 Kb and P < 5e-10 cutoff    7  rs62481385  rs67338333  1
    ## 136 10 Kb and P < 5e-10 cutoff    7   rs3757780  rs12668077  1
    ## 137 10 Kb and P < 5e-10 cutoff    8    rs268561    rs268624  1
    ## 138 10 Kb and P < 5e-10 cutoff    8    rs268561    rs268572  1
    ## 139 10 Kb and P < 5e-10 cutoff    8    rs268624    rs268572  1
    ## 140 10 Kb and P < 5e-10 cutoff    8    rs268561    rs187907  1
    ## 141 10 Kb and P < 5e-10 cutoff    8    rs268624    rs187907  1
    ## 142 10 Kb and P < 5e-10 cutoff    8    rs268572    rs187907  1
    ## 143 10 Kb and P < 5e-10 cutoff    8   rs7001162  rs11985375  1
    ## 144 10 Kb and P < 5e-10 cutoff    8  rs13280922  rs13281864  1
    ## 145 10 Kb and P < 5e-10 cutoff    8  rs13280922   rs7016334  1
    ## 146 10 Kb and P < 5e-10 cutoff    8  rs13281864   rs7016334  1
    ## 147 10 Kb and P < 5e-10 cutoff    8  rs13280922   rs6991708  1
    ## 148 10 Kb and P < 5e-10 cutoff    8  rs13281864   rs6991708  1
    ## 149 10 Kb and P < 5e-10 cutoff    8   rs7016334   rs6991708  1
    ## 150 10 Kb and P < 5e-10 cutoff    8  rs13280922   rs1838392  1
    ## 151 10 Kb and P < 5e-10 cutoff    8  rs13281864   rs1838392  1
    ## 152 10 Kb and P < 5e-10 cutoff    8   rs7016334   rs1838392  1
    ## 153 10 Kb and P < 5e-10 cutoff    8   rs6991708   rs1838392  1
    ## 154 10 Kb and P < 5e-10 cutoff    8  rs13280922  rs10429277  1
    ## 155 10 Kb and P < 5e-10 cutoff    8  rs13281864  rs10429277  1
    ## 156 10 Kb and P < 5e-10 cutoff    8   rs7016334  rs10429277  1
    ## 157 10 Kb and P < 5e-10 cutoff    8   rs6991708  rs10429277  1
    ## 158 10 Kb and P < 5e-10 cutoff    8   rs1838392  rs10429277  1
    ## 159 10 Kb and P < 5e-10 cutoff    8  rs13280922  rs11994908  1
    ## 160 10 Kb and P < 5e-10 cutoff    8  rs13281864  rs11994908  1
    ## 161 10 Kb and P < 5e-10 cutoff    8   rs7016334  rs11994908  1
    ## 162 10 Kb and P < 5e-10 cutoff    8   rs6991708  rs11994908  1
    ## 163 10 Kb and P < 5e-10 cutoff    8   rs1838392  rs11994908  1
    ## 164 10 Kb and P < 5e-10 cutoff    8  rs10429277  rs11994908  1
    ## 165 10 Kb and P < 5e-10 cutoff    8  rs13280922  rs13273979  1
    ## 166 10 Kb and P < 5e-10 cutoff    8  rs13281864  rs13273979  1
    ## 167 10 Kb and P < 5e-10 cutoff    8   rs7016334  rs13273979  1
    ## 168 10 Kb and P < 5e-10 cutoff    8   rs6991708  rs13273979  1
    ## 169 10 Kb and P < 5e-10 cutoff    8   rs1838392  rs13273979  1
    ## 170 10 Kb and P < 5e-10 cutoff    8  rs10429277  rs13273979  1
    ## 171 10 Kb and P < 5e-10 cutoff    8  rs11994908  rs13273979  1
    ## 172 10 Kb and P < 5e-10 cutoff    8  rs13280922  rs13255849  1
    ## 173 10 Kb and P < 5e-10 cutoff    8  rs13281864  rs13255849  1
    ## 174 10 Kb and P < 5e-10 cutoff    8   rs7016334  rs13255849  1
    ## 175 10 Kb and P < 5e-10 cutoff    8   rs6991708  rs13255849  1
    ## 176 10 Kb and P < 5e-10 cutoff    8   rs1838392  rs13255849  1
    ## 177 10 Kb and P < 5e-10 cutoff    8  rs10429277  rs13255849  1
    ## 178 10 Kb and P < 5e-10 cutoff    8  rs11994908  rs13255849  1
    ## 179 10 Kb and P < 5e-10 cutoff    8  rs13273979  rs13255849  1
    ## 180 10 Kb and P < 5e-10 cutoff    8   rs7001162  rs13274381  1
    ## 181 10 Kb and P < 5e-10 cutoff    8  rs11985375  rs13274381  1
    ## 182 10 Kb and P < 5e-10 cutoff    8  rs13280922   rs4446768  1
    ## 183 10 Kb and P < 5e-10 cutoff    8  rs13281864   rs4446768  1
    ## 184 10 Kb and P < 5e-10 cutoff    8   rs7016334   rs4446768  1
    ## 185 10 Kb and P < 5e-10 cutoff    8   rs6991708   rs4446768  1
    ## 186 10 Kb and P < 5e-10 cutoff    8   rs1838392   rs4446768  1
    ## 187 10 Kb and P < 5e-10 cutoff    8  rs10429277   rs4446768  1
    ## 188 10 Kb and P < 5e-10 cutoff    8  rs11994908   rs4446768  1
    ## 189 10 Kb and P < 5e-10 cutoff    8  rs13273979   rs4446768  1
    ## 190 10 Kb and P < 5e-10 cutoff    8  rs13255849   rs4446768  1
    ## 191 10 Kb and P < 5e-10 cutoff    8  rs13280922   rs7003794  1
    ## 192 10 Kb and P < 5e-10 cutoff    8  rs13281864   rs7003794  1
    ## 193 10 Kb and P < 5e-10 cutoff    8   rs7016334   rs7003794  1
    ## 194 10 Kb and P < 5e-10 cutoff    8   rs6991708   rs7003794  1
    ## 195 10 Kb and P < 5e-10 cutoff    8   rs1838392   rs7003794  1
    ## 196 10 Kb and P < 5e-10 cutoff    8  rs10429277   rs7003794  1
    ## 197 10 Kb and P < 5e-10 cutoff    8  rs11994908   rs7003794  1
    ## 198 10 Kb and P < 5e-10 cutoff    8  rs13273979   rs7003794  1
    ## 199 10 Kb and P < 5e-10 cutoff    8  rs13255849   rs7003794  1
    ## 200 10 Kb and P < 5e-10 cutoff    8   rs4446768   rs7003794  1
    ## 201 10 Kb and P < 5e-10 cutoff    8  rs13280922  rs13255749  1
    ## 202 10 Kb and P < 5e-10 cutoff    8  rs13281864  rs13255749  1
    ## 203 10 Kb and P < 5e-10 cutoff    8   rs7016334  rs13255749  1
    ## 204 10 Kb and P < 5e-10 cutoff    8   rs6991708  rs13255749  1
    ## 205 10 Kb and P < 5e-10 cutoff    8   rs1838392  rs13255749  1
    ## 206 10 Kb and P < 5e-10 cutoff    8  rs10429277  rs13255749  1
    ## 207 10 Kb and P < 5e-10 cutoff    8  rs11994908  rs13255749  1
    ## 208 10 Kb and P < 5e-10 cutoff    8  rs13273979  rs13255749  1
    ## 209 10 Kb and P < 5e-10 cutoff    8  rs13255849  rs13255749  1
    ## 210 10 Kb and P < 5e-10 cutoff    8   rs4446768  rs13255749  1
    ## 211 10 Kb and P < 5e-10 cutoff    8   rs7003794  rs13255749  1
    ## 212 10 Kb and P < 5e-10 cutoff    8  rs13280922   rs4581086  1
    ## 213 10 Kb and P < 5e-10 cutoff    8  rs13281864   rs4581086  1
    ## 214 10 Kb and P < 5e-10 cutoff    8   rs7016334   rs4581086  1
    ## 215 10 Kb and P < 5e-10 cutoff    8   rs6991708   rs4581086  1
    ## 216 10 Kb and P < 5e-10 cutoff    8   rs1838392   rs4581086  1
    ## 217 10 Kb and P < 5e-10 cutoff    8  rs10429277   rs4581086  1
    ## 218 10 Kb and P < 5e-10 cutoff    8  rs11994908   rs4581086  1
    ## 219 10 Kb and P < 5e-10 cutoff    8  rs13273979   rs4581086  1
    ## 220 10 Kb and P < 5e-10 cutoff    8  rs13255849   rs4581086  1
    ## 221 10 Kb and P < 5e-10 cutoff    8   rs4446768   rs4581086  1
    ## 222 10 Kb and P < 5e-10 cutoff    8   rs7003794   rs4581086  1
    ## 223 10 Kb and P < 5e-10 cutoff    8  rs13255749   rs4581086  1
    ## 224 10 Kb and P < 5e-10 cutoff    8  rs13280922   rs6984663  1
    ## 225 10 Kb and P < 5e-10 cutoff    8  rs13281864   rs6984663  1
    ## 226 10 Kb and P < 5e-10 cutoff    8   rs7016334   rs6984663  1
    ## 227 10 Kb and P < 5e-10 cutoff    8   rs6991708   rs6984663  1
    ## 228 10 Kb and P < 5e-10 cutoff    8   rs1838392   rs6984663  1
    ## 229 10 Kb and P < 5e-10 cutoff    8  rs10429277   rs6984663  1
    ## 230 10 Kb and P < 5e-10 cutoff    8  rs11994908   rs6984663  1
    ## 231 10 Kb and P < 5e-10 cutoff    8  rs13273979   rs6984663  1
    ## 232 10 Kb and P < 5e-10 cutoff    8  rs13255849   rs6984663  1
    ## 233 10 Kb and P < 5e-10 cutoff    8   rs4446768   rs6984663  1
    ## 234 10 Kb and P < 5e-10 cutoff    8   rs7003794   rs6984663  1
    ## 235 10 Kb and P < 5e-10 cutoff    8  rs13255749   rs6984663  1
    ## 236 10 Kb and P < 5e-10 cutoff    8   rs4581086   rs6984663  1
    ## 237 10 Kb and P < 5e-10 cutoff    8   rs7001162  rs56289931  1
    ## 238 10 Kb and P < 5e-10 cutoff    8  rs11985375  rs56289931  1
    ## 239 10 Kb and P < 5e-10 cutoff    8  rs13274381  rs56289931  1
    ## 240 10 Kb and P < 5e-10 cutoff    8  rs13280922   rs2008517  1
    ## 241 10 Kb and P < 5e-10 cutoff    8  rs13281864   rs2008517  1
    ## 242 10 Kb and P < 5e-10 cutoff    8   rs7016334   rs2008517  1
    ## 243 10 Kb and P < 5e-10 cutoff    8   rs6991708   rs2008517  1
    ## 244 10 Kb and P < 5e-10 cutoff    8   rs1838392   rs2008517  1
    ## 245 10 Kb and P < 5e-10 cutoff    8  rs10429277   rs2008517  1
    ## 246 10 Kb and P < 5e-10 cutoff    8  rs11994908   rs2008517  1
    ## 247 10 Kb and P < 5e-10 cutoff    8  rs13273979   rs2008517  1
    ## 248 10 Kb and P < 5e-10 cutoff    8  rs13255849   rs2008517  1
    ## 249 10 Kb and P < 5e-10 cutoff    8   rs4446768   rs2008517  1
    ## 250 10 Kb and P < 5e-10 cutoff    8   rs7003794   rs2008517  1
    ## 251 10 Kb and P < 5e-10 cutoff    8  rs13255749   rs2008517  1
    ## 252 10 Kb and P < 5e-10 cutoff    8   rs4581086   rs2008517  1
    ## 253 10 Kb and P < 5e-10 cutoff    8   rs6984663   rs2008517  1
    ## 254 10 Kb and P < 5e-10 cutoff    8  rs13280922   rs2732143  1
    ## 255 10 Kb and P < 5e-10 cutoff    8  rs13281864   rs2732143  1
    ## 256 10 Kb and P < 5e-10 cutoff    8   rs7016334   rs2732143  1
    ## 257 10 Kb and P < 5e-10 cutoff    8   rs6991708   rs2732143  1
    ## 258 10 Kb and P < 5e-10 cutoff    8   rs1838392   rs2732143  1
    ## 259 10 Kb and P < 5e-10 cutoff    8  rs10429277   rs2732143  1
    ## 260 10 Kb and P < 5e-10 cutoff    8  rs11994908   rs2732143  1
    ## 261 10 Kb and P < 5e-10 cutoff    8  rs13273979   rs2732143  1
    ## 262 10 Kb and P < 5e-10 cutoff    8  rs13255849   rs2732143  1
    ## 263 10 Kb and P < 5e-10 cutoff    8   rs4446768   rs2732143  1
    ## 264 10 Kb and P < 5e-10 cutoff    8   rs7003794   rs2732143  1
    ## 265 10 Kb and P < 5e-10 cutoff    8  rs13255749   rs2732143  1
    ## 266 10 Kb and P < 5e-10 cutoff    8   rs4581086   rs2732143  1
    ## 267 10 Kb and P < 5e-10 cutoff    8   rs6984663   rs2732143  1
    ## 268 10 Kb and P < 5e-10 cutoff    8   rs2008517   rs2732143  1
    ## 269 10 Kb and P < 5e-10 cutoff    8  rs13280922   rs2732127  1
    ## 270 10 Kb and P < 5e-10 cutoff    8  rs13281864   rs2732127  1
    ## 271 10 Kb and P < 5e-10 cutoff    8   rs7016334   rs2732127  1
    ## 272 10 Kb and P < 5e-10 cutoff    8   rs6991708   rs2732127  1
    ## 273 10 Kb and P < 5e-10 cutoff    8   rs1838392   rs2732127  1
    ## 274 10 Kb and P < 5e-10 cutoff    8  rs10429277   rs2732127  1
    ## 275 10 Kb and P < 5e-10 cutoff    8  rs11994908   rs2732127  1
    ## 276 10 Kb and P < 5e-10 cutoff    8  rs13273979   rs2732127  1
    ## 277 10 Kb and P < 5e-10 cutoff    8  rs13255849   rs2732127  1
    ## 278 10 Kb and P < 5e-10 cutoff    8   rs4446768   rs2732127  1
    ## 279 10 Kb and P < 5e-10 cutoff    8   rs7003794   rs2732127  1
    ## 280 10 Kb and P < 5e-10 cutoff    8  rs13255749   rs2732127  1
    ## 281 10 Kb and P < 5e-10 cutoff    8   rs4581086   rs2732127  1
    ## 282 10 Kb and P < 5e-10 cutoff    8   rs6984663   rs2732127  1
    ## 283 10 Kb and P < 5e-10 cutoff    8   rs2008517   rs2732127  1
    ## 284 10 Kb and P < 5e-10 cutoff    8   rs2732143   rs2732127  1
    ## 285 10 Kb and P < 5e-10 cutoff    8  rs13280922   rs2639942  1
    ## 286 10 Kb and P < 5e-10 cutoff    8  rs13281864   rs2639942  1
    ## 287 10 Kb and P < 5e-10 cutoff    8   rs7016334   rs2639942  1
    ## 288 10 Kb and P < 5e-10 cutoff    8   rs6991708   rs2639942  1
    ## 289 10 Kb and P < 5e-10 cutoff    8   rs1838392   rs2639942  1
    ## 290 10 Kb and P < 5e-10 cutoff    8  rs10429277   rs2639942  1
    ## 291 10 Kb and P < 5e-10 cutoff    8  rs11994908   rs2639942  1
    ## 292 10 Kb and P < 5e-10 cutoff    8  rs13273979   rs2639942  1
    ## 293 10 Kb and P < 5e-10 cutoff    8  rs13255849   rs2639942  1
    ## 294 10 Kb and P < 5e-10 cutoff    8   rs4446768   rs2639942  1
    ## 295 10 Kb and P < 5e-10 cutoff    8   rs7003794   rs2639942  1
    ## 296 10 Kb and P < 5e-10 cutoff    8  rs13255749   rs2639942  1
    ## 297 10 Kb and P < 5e-10 cutoff    8   rs4581086   rs2639942  1
    ## 298 10 Kb and P < 5e-10 cutoff    8   rs6984663   rs2639942  1
    ## 299 10 Kb and P < 5e-10 cutoff    8   rs2008517   rs2639942  1
    ## 300 10 Kb and P < 5e-10 cutoff    8   rs2732143   rs2639942  1
    ## 301 10 Kb and P < 5e-10 cutoff    8   rs2732127   rs2639942  1
    ## 302 10 Kb and P < 5e-10 cutoff    8   rs3098869   rs3110262  1
    ## 303 10 Kb and P < 5e-10 cutoff    8   rs3098869   rs4074908  1
    ## 304 10 Kb and P < 5e-10 cutoff    8   rs3110262   rs4074908  1
    ## 305 10 Kb and P < 5e-10 cutoff    9  rs10793962  rs75444660  1
    ## 306 10 Kb and P < 5e-10 cutoff    9  rs78590974  rs28632066  1

Here we calculate average r<sup>2</sup> value among the SNPs (IVs).

``` r
ld_mat_tb %>% 
  group_by(dist_cutoff) %>% 
  summarise(avg_r2 = mean(r2, na.rm = TRUE))
```

    ## # A tibble: 2 × 2
    ##   dist_cutoff                    avg_r2
    ##   <fct>                           <dbl>
    ## 1 "10 Kb and P < 5e-10 cutoff " 0.274  
    ## 2 "1 Mb and P < 5e-8 cutoff"    0.00558

Here are the top 3 highest r<sup>2</sup> values, when distance cutoff
was increased to 1 Mb or 10 Mb.

``` r
ld_mat_tb %>% 
  filter(!grepl("Kb", dist_cutoff)) %>% 
  group_by(dist_cutoff) %>% 
  top_n(n=3, wt = r2)
```

    ## # A tibble: 3 × 5
    ## # Groups:   dist_cutoff [1]
    ##   dist_cutoff              chr   SNP1      SNP2          r2
    ##   <fct>                    <fct> <fct>     <fct>      <dbl>
    ## 1 1 Mb and P < 5e-8 cutoff 2     rs4671051 rs7423637  0.038
    ## 2 1 Mb and P < 5e-8 cutoff 2     rs4671051 rs13017210 0.056
    ## 3 1 Mb and P < 5e-8 cutoff 7     rs7795564 rs3757582  0.049

 

In conclusion, by using the thresholds described by Zhu *et al.*, we
show that multiple genetic variants are in strong LD, and thus
correlating (that is, are not independent). Using higher distance
cutoffs

 

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
    ##  [1] LDlinkR_1.3.0     TwoSampleMR_0.5.6 patchwork_1.1.2   forcats_0.5.2    
    ##  [5] stringr_1.4.1     dplyr_1.1.2       purrr_0.3.4       readr_2.1.2      
    ##  [9] tidyr_1.2.1       tibble_3.2.1      ggplot2_3.4.2     tidyverse_1.3.2  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] httr_1.4.4          viridisLite_0.4.1   jsonlite_1.8.0     
    ##  [4] splines_4.2.1       foreach_1.5.2       R.utils_2.12.2     
    ##  [7] modelr_0.1.9        assertthat_0.2.1    highr_0.9          
    ## [10] googlesheets4_1.0.1 cellranger_1.1.0    yaml_2.3.5         
    ## [13] pillar_1.9.0        backports_1.4.1     lattice_0.20-45    
    ## [16] glue_1.6.2          digest_0.6.29       rvest_1.0.3        
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
    ## [64] DBI_1.1.3           curl_4.3.2          reshape2_1.4.4     
    ## [67] R6_2.5.1            lubridate_1.9.2     knitr_1.40         
    ## [70] fastmap_1.1.0       utf8_1.2.2          nortest_1.0-4      
    ## [73] shape_1.4.6         stringi_1.7.8       Rcpp_1.0.9         
    ## [76] vctrs_0.6.3         dbplyr_2.2.1        tidyselect_1.2.0   
    ## [79] xfun_0.33
