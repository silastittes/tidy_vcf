# tidy_vcf
Converts a sample of sites from standard VCF table into tidy format for easy summary stats and plotting.

It's often desirable to get a sense of the site and genotype metrics from a VCF file prior to filtering. This simple program is designed to faciltate that. While there are excellent tools like [`bcftools stats`](http://samtools.github.io/bcftools/bcftools.html#stats) to do this already, the summary statistics gathered are very simple and not very customizable. If you don't mind exchanging speed for more information and flexibility, give this tool a try.

The program determines the data frame features based on the FORMAT and INFO listed in the VCF header. The output features depend on what is listed in the header. When an individual site does not contain a particular INFO feature, an `NA` is added to the row of that feature. This means the output data frames should contain all the information in the original. It also means the VCF header that faithfully represents what is present in the VCF rows is *absolutely essential.*

Depends on `python3` with the `argparse`, `gzip`, and `random` libraries available.  


## Installation

`pip install tidy-vcf`

## Usage

```
usage: tidy_vcf [-h] [-s [SITES]] [-t THIN] -v [VCF] -o [SITES_OUT] -g
                [GENOTYPE_OUT]

Given a vcf file, produces a tidy versions of sites and genotypes data, in
efforts to make it easier to calculate summary statitics and visualizations.

optional arguments:
  -h, --help            show this help message and exit
  -s [SITES], --sites [SITES]
                        Input file listing tab delimited sites to output. Each
                        line must be a chromosome/scaffold name and the
                        1-indexed position of the site.
  -t THIN, --thin THIN  Alternative to --sites, where a sites will be selected
                        no less than THIN bases apart.
  -v [VCF], --vcf [VCF]
                        VCF input file. Should end in .vcf or .gz (plain text
                        or compressed).
  -o [SITES_OUT], --sites_out [SITES_OUT]
                        Name of output file for tidy sites data. Cannot be
                        used simultaneously with --thin
  -g [GENOTYPE_OUT], --genotype_out [GENOTYPE_OUT]
                        Name of output file for tidy genotype data.
```


## Working example command:

This repo contains necessary examples to test out the utility, but three example commands assume you have installed the program using `pip` above. 

For example `sites.txt` contains 1,000 sites, the following command will generate to two tidy data files:

```
tidy_vcf -s example/sites.txt -v example/test.vcf.gz -o example/test_tidy_sites.vcf -g example/test_tidy_genotypes.vcf
```

Alternatively, you can also using the `-t` argument, which specifies how far apart sampled positions should be:

```
tidy_vcf -t 1000 -v example/test.vcf.gz -o example/test_tidy_sites.vcf -g example/test_tidy_genotypes.vcf
```


Alternatively, you can use neither `-t` or `-s`, in which case all sites will be used with a warning: 

```
tidy_vcf -v example/test.vcf.gz -o example/test_tidy_sites.vcf -g example/test_tidy_genotypes.vcf
```


# Now what?

The sites and genotypes file are easy to load into R or Python for summary and visualzation.

Here's a simple example that can be done in R:

```
#Read sites data
sites <- vroom::vroom("example/test_tidy_sites.vcf", delim = "\t", na = c("NA", "./.", '.'), comment = "#")

head(sites)
```

```
# A tibble: 6 x 23
  CHROM   POS ID    REF   ALT     QUAL FILTER    AC AF       AN BaseQRankSum ClippingRankSum    DP ExcessHet
  <chr> <dbl> <lgl> <chr> <chr>  <dbl> <lgl>  <dbl> <chr> <dbl>        <dbl>           <dbl> <dbl>     <dbl>
1 chr1   4240 NA    TC    T,CC   302.  NA       113 0.01…    62        0.277               0   179    3.13  
2 chr1   4248 NA    A     G     3303.  NA        68 1.000    68       NA                  NA   186    3.01  
3 chr1   4261 NA    C     T      113.  NA        15 0.214    70       -3.37                0   218    0     
4 chr1   4375 NA    C     T       41.3 NA         2 0.007   284        0.089               0   878    0.0077
5 chr1   4441 NA    G     A      110.  NA         2 0.021    96        0.713               0   212    0.0229
6 chr1   4495 NA    C     T      280.  NA        21 0.250    84        2.07                0   130    0.0017
# … with 9 more variables: FS <dbl>, InbreedingCoeff <dbl>, MLEAC <dbl>, MLEAF <chr>, MQ <dbl>, MQRankSum <dbl>,
#   QD <dbl>, ReadPosRankSum <dbl>, SOR <dbl>
```

```
#Read genotype data
genos <-vroom::vroom("example/test_tidy_genotypes.vcf", delim = "\t", na = c("NA", "./.", '.'), comment = "#")  
head(genos)
```

```
# A tibble: 6 x 13
  IND        CHROM   POS ID    REF   ALT    QUAL FILTER GT    AD       DP    GQ PL   
  <chr>      <chr> <dbl> <lgl> <chr> <chr> <dbl> <lgl>  <chr> <chr> <dbl> <dbl> <chr>
1 JRIAL10-10 chr1   4240 NA    TC    T,CC   302. NA     NA    NA        0     0 NA   
2 JRIAL10-11 chr1   4240 NA    TC    T,CC   302. NA     NA    NA        0     0 NA   
3 JRIAL10-12 chr1   4240 NA    TC    T,CC   302. NA     NA    NA        0     0 NA   
4 JRIAL10-13 chr1   4240 NA    TC    T,CC   302. NA     NA    NA        0     0 NA   
5 JRIAL10-14 chr1   4240 NA    TC    T,CC   302. NA     NA    NA        0     0 NA   
6 JRIAL10-15 chr1   4240 NA    TC    T,CC   302. NA     NA    NA        0     0 NA 
```


```
#For each individual, summarize average proportion of missing genotypes, and average and standard deviation of  depth and genotype quality.
genos %>% 
  group_by(IND) %>% 
  summarise(
    missing = mean(is.na(GT)),
    dp = mean(DP),
    dp_std = sd(DP),
    gq = mean(GQ),
    gq_std = sd(GQ)
    )
```

```
# A tibble: 195 x 6
   IND       missing    dp dp_std    gq gq_std
   <chr>       <dbl> <dbl>  <dbl> <dbl>  <dbl>
 1 JRIAL1-1    0.706 1.05    2.44  3.71   9.80
 2 JRIAL1-10   0.753 0.662   1.69  2.37   7.03
 3 JRIAL1-11   0.633 1.77    3.93  5.67  13.0 
 4 JRIAL1-12   0.701 1.18    2.85  3.99  10.5 
 5 JRIAL1-13   0.662 1.31    2.92  4.49  11.1 
 6 JRIAL1-14   0.693 1.07    2.43  3.78   9.99
 7 JRIAL1-15   0.595 1.15    2.09  3.91   7.93
 8 JRIAL1-16   0.722 0.788   1.88  2.70   7.74
 9 JRIAL1-17   0.675 1.13    2.44  3.96   9.81
10 JRIAL1-18   0.714 0.848   2.05  2.96   8.00
# … with 185 more rows
```



## How to get sites

Unless a set of sites is of interest, it probably makes more sense to use the `--thin` flag. However, if you don't want to use `--thin` one way to get a random set would be:

```
zcat < example/test.vcf.gz | grep -v "#" | cut -f1,2 | shuf -n 100 | sort -V > sites.txt
```

where `shuf -n` could be modified to a desired number. If your VCF is very large, another approach to get a set would likely be needed. 

