# tidy_vcf
Converts a sample of sites from standard VCF table into tidy format for easy summary stats and plotting.

It's often desirable to get a sense of the site and genotype metrics from a VCF file prior to filtering. This simple program is designed to faciltate that. While there are excellent tools like `bcftools stats` to do this already, the summary statistics gathered are very simple and not very customizable. If you don't mind exchanging speed for more information and flexibility, give this tool a try.

```
usage: tidy_vcf.py [-h] -s [SITES] -v [VCF] -o [SITES_OUT] -g [GENOTYPE_OUT]

Given a vcf file and a file of sites, produces a tidy versions of sites and
genotypes data, in efforts to make it easier to calculate summary statitics
and visualizations.

optional arguments:
  -h, --help            show this help message and exit
  -s [SITES], --sites [SITES]
                        Input file listing tab delimited sites to output. Each
                        line must be a chromosome/scaffold name and the
                        1-indexed position of the site.
  -v [VCF], --vcf [VCF]
                        VCF input file. Should end in .vcf or .gz (plain text
                        or compressed).
  -o [SITES_OUT], --sites_out [SITES_OUT]
                        Name of output file for tidy sites data.
  -g [GENOTYPE_OUT], --genotype_out [GENOTYPE_OUT]
                        Name of output file for tidy genotype data.
```


#Working example command:

For example `sites.txt` contains 1,000 sites, the following command will generate to two tidy data files:

```
python tidy_vcf.py -s sites.txt -v test.vcf.gz -o test_tidy_sites.vcf -g test_tidy_genotypes.vcf
```

#Now what?

The sites and genotypes file are easy to load into R or Python for summary and visualzation.

Here's a simple example that can be done in R:

```
#Read sites data
sites <- vroom::vroom("test_tidy_sites.vcf", delim = "\t", na = c("NA", "./.", '.'), comment = "#")

#Read genotype data
genos <-vroom::vroom("test_tidy_genotypes.vcf", delim = "\t", na = c("NA", "./.", '.'), comment = "#")  

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



#How to get sites

One tedious part of this is needing sites to sample before hand. Assuming you don't already have interest in some set, one way to get a random set would be:

```
zcat < test.vcf.gz | grep -v "#" | cut -f1,2 | shuf -n 100 | sort -V > sites.txt
```

where `shuf -n` could be modified to a desired number. If your VCF is very large, another approach would likely be needed. 

