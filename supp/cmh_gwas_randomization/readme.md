### Find the number of overlaps between the CMH scan with FDR < 0.01 and the 893 (non clumped) GWAS loci

```
gwasbed=/scratch/midway2/rozennpineau/drought/compare_sites_commongarden_drought/drought/gwas_893.bed
cmhbed=/scratch/midway2/rozennpineau/drought/compare_sites_commongarden_drought/drought/CMH_FDR01.bed
bedtools intersect -a $gwasbed -b $cmhbed  -wa -wb -f 0.99 -r > overlap_cmh_gwas.txt
```
### Find the number of genes in the clumped and non-clumped CMH scan with FDR < 0.01

How many genes in the non-clumped CMH scan (FDR < 0.01) ?
```
cmhbed=/scratch/midway2/rozennpineau/drought/compare_sites_commongarden_drought/drought/CMH_FDR01.bed
#prep bed file
awk '{OFS="\t";print $1, $2, $2}' $cmhbed > cmh_010.bed

#add “Scaffold_” 
awk -v OFS="\t" '{$1 = "Scaffold_" $1}1' cmh_010.bed > cmh_010_scaffold.bed

#intersect with gff
act-conda
gff=/project/kreiner/data/genome/Atub_193_hap2.all.sorted.gff
bed=cmh_010_scaffold.bed
bedtools intersect -a $bed -b $gff -wo > overlap_cmh_010_gff.bed

#extract number of CMH hits in genes
less overlap_cmh_010_gff.bed | cut -f 6 | grep gene | wc -l
#44,692
```
44,692 hits in genes !


Clumpimpg the CMH file needs several steps:
**(1)** Filter initial vcf based on FDR < 0.01 bed file
```
#!/bin/bash
#SBATCH --job-name=filter_vcf
#SBATCH --output=err.out
#SBATCH --error=err.in
#SBATCH --time=5:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10G   # memory per cpu-core

cd /scratch/midway2/rozennpineau/drought/compare_sites_commongarden_drought/drought
module load vcftools
vcf=/scratch/midway3/rozennpineau/drought/admixture/merged_numericChr.vcf.gz
bed=/scratch/midway2/rozennpineau/drought/compare_sites_commongarden_drought/drought/CMH_FDR01.bed
vcftools --gzvcf $vcf --bed $bed --out cmh_010 --recode

#note : check if you need to use the character (Scaffold_) and not numeric names!
#note : need for a slurm job, the terminal without a batch job was very slow!
```
**(2)** Make plink files

```
#change CHROM:POS field to ID
module load htslib
bgzip cmh_010.recode.vcf
bcftools tabix cmh_010.recode.vcf.gz
bcftools annotate --set-id +'%CHROM:%POS' cmh_010.recode.vcf.gz -o cmh_010_ID.vcf.gz

#make plink family files
vcf=cmh_010_ID.vcf.gz
plink --vcf $vcf --out cmh_010_tothin --allow-extra-chr --recode --double-id
```

```
#in: /scratch/midway2/rozennpineau/drought/compare_sites_commongarden_drought/drought
assoc=FDRdrought 
#remove Scaffold from chromosome name in association file
sed 's/Scaffold_//g' $assoc > FDRdrought_numeric
#replace SNP with ID in association file
sed 's/SNP/ID/g' FDRdrought_numeric > FDRdrought_numeric_ID

assoc=FDRdrought_numeric_ID
plink --file cmh_010_tothin --clump $assoc --clump-p1 0.01 --clump-field FDR_p --clump-kb 100 --out cmh_010_clumped_100kb --allow-no-sex --allow-extra-chr --clump-snp-field ID

```
42831 clumps formed from 94052 top variants.

How many genes in the clumped CMH scan (FDR < 0.01) ?

```
cmhbed=/scratch/midway2/rozennpineau/drought/compare_sites_commongarden_drought/drought/CMH_FDR01.bed
#prep bed file
awk '{OFS="\t";print $1, $2, $2}' $cmhbed > cmh_010.bed

#add “Scaffold_” 
awk -v OFS="\t" '{$1 = "Scaffold_" $1}1' cmh_010.bed > cmh_010_scaffold.bed

#intersect with gff
act-conda
gff=/project/kreiner/data/genome/Atub_193_hap2.all.sorted.gff
bed=cmh_010_scaffold.bed
bedtools intersect -a $bed -b $gff -wo > overlap_cmh_010_gff.bed

#extract number of CMH hits in genes
less overlap_cmh_010_gff.bed | cut -f 6 | grep gene | wc -l
```



### Randomization analysis for CMH scan / drought variants taking into account LD

(1) extract bed file that spans mean ancestry tract length for each 43 lead SNP, with the lead SNP placed in the middle of the tract:

```
#Extract mean length ancestry block for each 43 lead SNPs
lead_snp_bed <- read.table("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/8.ancestry_hmm/1.results/drought_adapted_loci_43.txt", sep = "\t", header = T)
lead_snp_bed <- lead_snp_bed[,c(1,2)]
lead_snp_bed$len <- NA

#upload the mean ancestry tract length per 
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/11.common_garden/randomization/")
mean_len_genome <- read.table("full_genome_anc_block_mean_len_tab.bed", sep = "\t", header = T)

for (snp in 1:43) {
  
  idx <- which(lead_snp_bed$chrom[snp] == mean_len_genome$chrom & lead_snp_bed$pos[snp] == mean_len_genome$pos1)
  lead_snp_bed$len <- mean_len_genome$mean_anc_block_length[idx]
  
}

#Make bed file with mean length spanning those 43 lead snps
lead_snp_mean_clumps <- data.frame(chrom = lead_snp_bed$chrom )
lead_snp_mean_clumps$pos1 <- lead_snp_bed$pos - round(lead_snp_bed$len/2)
lead_snp_mean_clumps$pos2 <- lead_snp_bed$pos + round(lead_snp_bed$len/2)

#Export for intersection with CMH
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/11.common_garden/randomization/")
write.table(lead_snp_mean_clumps, "lead_snp_mean_clumps.bed", sep = "\t", col.names = T, row.names = F, quote = F)
```
**Note : do we want to keep SNPs that are at the edge between two ancestry tracts ?
--> maybe I need to re-work the script that calculates the lengths of each tract to make sure it starts at index + 1 for the next tract

(2) count the number of CMH SNPs that fall within those ancestry tracts.

```
bash

gwasbed=/scratch/midway3/rozennpineau/drought/randomization/1_overlap_cmh_gwas_observed/lead_snp_mean_clumps.bed 
cmhbed=/scratch/midway2/rozennpineau/drought/compare_sites_commongarden_drought/drought/CMH_FDR01.bed

bedtools intersect -a $gwasbed -b $cmhbed -c > intersect_cmh_43meanblock.txt
#-c For each feature in A, report the number of overlapping features in B

#record the total and mean number of overlaps
awk '{ sum += $NF; n++ } END { print "Total:", sum, "Mean:", sum/n }' intersect_cmh_43meanblock.txt
#Total: 8381 Mean: 194.907

```
I have the total number of SNPs from the CMH that fall within the ancestry tracts from the 43 lead SNPs. 
Next is to randomly sample bed file that have the same mean size distribution of clumps, and record the mean and total number for those.

(3) create the randomly sampled bed files:

```
r

#goal : to create bed files that have a similar distribution of ancestry block lengths as the observed drought bed file

#full <- read.table("../full_genome_anc_block_mean_len_tab_subset.bed", sep = "\t", header = F) #change to header = T with full file
full <- read.table("../ancestry_block_mean_786257.bed", sep = "\t", header = T)
#obs <- read.table("observed_anc_block_len_distr.txt", sep = "\t", header = T) #plus or minus 100 bp
obs <- read.table("../1_overlap_cmh_gwas_observed/lead_snp_mean_clumps.bed", sep = "\t", header = T ) #43 clumps
n_perm <- 100
precision <- 500
#precision <- 1000
precision <- 5000

for (perm in 1:n_perm) {
        cur_bed <- data.frame(matrix(NA, length(obs$len), 4))
        colnames(cur_bed) <- c("chrom", "pos1", "pos2", "len")
        for (i in 1:length(obs$len)) {
                idx_pool <- which(full[,4] > obs$len[i]-precision & full[,4] < obs$len[i]+precision)
                if (length(idx_pool)<1){break} #no matches
                idx_local <- sample(idx_pool,1,replace=F)
                cur_bed$chrom[i] <- full$chrom[idx_local]
                cur_bed$pos1[i] <- full$pos1[idx_local] - round(full$mean_anc_block_length[idx_local]/2)
                cur_bed$pos2[i] <- full$pos1[idx_local] + round(full$mean_anc_block_length[idx_local]/2)
                cur_bed$len[i] <- full$mean_anc_block_length[idx_local]
                write.table(cur_bed, paste("random_anc_block_", perm,".bed", sep = ""), sep = "\t", col.names = F, row.names = F, quote = F )
        }
}

```

(4) calculate the overlaps

```
bash


#!/bin/bash/
#SBATCH --job-name=rdn
#SBATCH --output=rdn.out
#SBATCH --error=rdn.err
#SBATCH --time=1:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10G   # memory per cpu-core

#activate conda
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner

cd /scratch/midway3/rozennpineau/drought/randomization/3_calculate_overlap/

#load cmh bed
cmhbed=/scratch/midway2/rozennpineau/drought/compare_sites_commongarden_drought/drought/CMH_FDR01.bed

rm -f nb_overlaps.txt
echo -e "total_overlaps\tmean_overlaps" > nb_overlaps.txt

for gwasbed in ../2_create_bed/*.bed; do

    bedtools intersect -a $gwasbed -b $cmhbed -c |
    awk '
        { sum += $NF; n++ }
        END {
            if (n > 0)
                printf "%d\t%.6f\n", sum, sum/n
            else
                printf "0\tNA\n"
        }
    ' >> nb_overlaps.txt

done
```


Plot QQplot and Manhattan plots: [plot_qqplot_cmh.R](https://github.com/rozenn-pineau/drought-project/edit/main/supp/cmh_gwas_randomization/plot_qqplot_cmh.R)

Plot randomization results: [plot_CMH_gwas.Rmd](https://github.com/rozenn-pineau/drought-project/edit/main/supp/cmh_gwas_randomization/plot_CMH_gwas.Rmd)















