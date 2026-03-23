### Calculating Tajima's D for the 43 loci under fluctuating selection 

```
# 0.5kb windows
vcftools --gzvcf drought_adapted_43clumps.recode.vcf.gz --TajimaD 10000 --remove-indv P16_Nat_1_T --remove-indv P12_Nat_14_T --out drought_adapted_43clumps_500bpwindows


# 10 kb windows
vcftools --gzvcf drought_adapted_43clumps.recode.vcf.gz --TajimaD 10000 --remove-indv P16_Nat_1_T --remove-indv P12_Nat_14_T --out drought_adapted_43clumps_10kbwindows


# 100 kb windows
vcftools --gzvcf drought_adapted_43clumps.recode.vcf.gz --TajimaD 100000 --remove-indv P16_Nat_1_T --remove-indv P12_Nat_14_T --out drought_adapted_43clumps_100kbwindows

```
### Compare to genome-wide signal 

```
# 0.5kb windows
vcftools --vcf two_pulse_flexible_prop_2_values.vcf --remove-indv P16_Nat_1_T --remove-indv P12_Nat_14_T --TajimaD 500 --out genome-wide_500bpwindows

# 10 kb windows, remove the two outliers from vcf
vcftools --vcf two_pulse_flexible_prop_2_values.vcf --remove-indv P16_Nat_1_T --remove-indv P12_Nat_14_T --TajimaD 10000 --out genome-wide_10kbwindows

# 100 kb windows, remove the two outliers from vcf
vcftools --vcf two_pulse_flexible_prop_2_values.vcf --remove-indv P16_Nat_1_T --remove-indv P12_Nat_14_T --TajimaD 100000 --out genome-wide_100kbwindows

```

Calculating Tajima's D on ancestry actually might not work (ancestry choice is 1 or 0, when choice for SNP can be 1 out of 4 different choices). Going back to the genotype file to calculate ancestry on the genotype:
### Filter vcf from drought bed (43 loci)

```
```

### Compare to genome-wide signal 

```
# 0.5kb windows
vcftools --gzvcf ../merged_filteredsnps_rename.vcf.gz --remove-indv P16_Nat_1_T --remove-indv P12_Nat_14_T --TajimaD 500 --out genome-wide_500bpwindows

# 10 kb windows, remove the two outliers from vcf
vcftools --gzvcf ../merged_filteredsnps_rename.vcf.gz --remove-indv P16_Nat_1_T --remove-indv P12_Nat_14_T --TajimaD 10000 --out genome-wide_10kbwindows

```


Plot distributions in R. 

