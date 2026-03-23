### Calculate Tajima's D for the 43 loci under fluctuating selection 

```
# 10 kb windows
vcftools --gzvcf drought_adapted_43clumps.recode.vcf.gz --TajimaD 10000 --out drought_adapted_43clumps_10kbwindows


# 100 kb windows
vcftools --gzvcf drought_adapted_43clumps.recode.vcf.gz --TajimaD 100000 --out drought_adapted_43clumps_100kbwindows

```
