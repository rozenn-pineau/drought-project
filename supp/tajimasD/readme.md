### Calculating Tajima's D for the 43 loci under fluctuating selection 

```
vcf=/scratch/midway2/rozennpineau/drought/two_pulse_flexible_prop_2/drought_adapted_43clumps.recode.vcf.gz

# 5 kb windows
cd /scratch/midway3/rozennpineau/drought/tajimasD/drought-adapted_loci/
vcftools --gzvcf $vcf --TajimaD 5000 --remove-indv P16_Nat_1_T --remove-indv P12_Nat_14_T --out drought_adapted_43clumps_5kbwindows


```
### Compare to genome-wide signal - randomization

(1) generate 43 loci long random bed files from the ancestry file
```
#goal : to create bed file by randomly subsmapling 43 loci across the genome

n_perm <- 1000
n_loci <- 43
bed <- read.table("drought_adapted_43clumps.bed", sep = "\t", header = T)
ancestry <- read.table("/cds3/kreiner/drought/analyses/randomization/ancestry_block_mean_786257.bed", sep = "\t", header = T)
        
for (perm in 1:n_perm) {
        cur_bed <- ancestry[sample(1:dim(ancestry)[1], n_loci),1:3]
        write.table(cur_bed, paste("random_43_loci", perm,".bed", sep = ""), sep = "\t", col.names = F, row.names = F, quote = F )
} 
```

(2) filter vcf from bed

```
vcf=/scratch/midway2/rozennpineau/drought/two_pulse_flexible_prop_2/two_pulse_flexible_prop_2_values.vcf
bedtools intersect -b random_43_loci1.bed -a $vcf > tmp1
#add header
cat header tmp1 > tmp2
vcftools --vcf tmp2 --TajimaD 5000 --remove-indv P16_Nat_1_T --remove-indv P12_Nat_14_T --out drought_adapted_43clumps_5kbwindows_test

```


