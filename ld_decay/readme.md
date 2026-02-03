### Estimating LD decay rate

We used PopLDdecay to estimate the LD decay rate on the genotypes (https://github.com/BGI-shenzhen/PopLDdecay). 
To decrease computational time, we initially subsampled the file by filtering for MAF <0.05, thinning one SNP every 2kb regardless of LD and keeping the variants only if present in >80% of samples.

```

#Run of MAF < 0.05 file
# Thin to keep 1 SNP every 2kb regardless of LD to decrease run time
plink \
  --bfile merged_numericChr_maf05 \
  --indep-pairwise 2000 1 0.99 \ #no LD filtering
  --out merged_numericChr_maf05_thin2kb
#running in a batch job because it was taking too long

#extract thinned SNPs from initial file
plink \
  --bfile merged_numericChr_maf05 \
  --extract merged_numericChr_maf05_thin2kb.prune.in \
  --make-bed \
  --out merged_numericChr_maf05_thin2kb

#Keep only if present in >80% of samples
plink \
  --bfile merged_numericChr_maf05_thin2kb \
  --geno 0.2 \
  --make-bed \
  --out merged_numericChr_maf05_thin2kb_miss80

#convert back to vcf for PopLDdecay
plink \
  --bfile merged_numericChr_maf05_thin2kb_miss80 \
  --recode vcf bgz \
  --out merged_numericChr_maf05_thin2kb_miss80


#run on vcf file directly
cd /scratch/midway2/rozennpineau/drought/compare_sites_commongarden_drought/drought/ld_decay/PopLDdecay/
./bin/PopLDdecay -InVCF /scratch/midway3/rozennpineau/drought/admixture/merged_numericChr_maf05_thin2kb_miss80.vcf.gz -OutStat LDdecay 
#figure for one population
perl  bin/Plot_OnePop.pl  -inFile   LDdecay.stat.gz  -output  Fig


```
