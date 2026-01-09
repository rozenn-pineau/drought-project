### CMH scan - GWAS output comparison
we want to compare the outputs from the CMH scans that identidied "agriculturally-adapted alleles" to the output of the GWAS on drought. I do this by comparing bed files both both (1) all loci, and (2) significant loci in both. 


On the cluster, I prepare the bed files for the CMH scans for comparison. 
```
#prepare the full CMH file bed (FDRdrought is the CMH scan output)
cat FDRdrought | \sed s/^\Scaffold_//g | awk -v OFS="\t" '{ print $1,$3,$3,$13}' > CMH_all.bed #63,979,747
#remove header
tail -n +2 CMH_all.bed > CMH_49338567.bed

#prepare the significant loci bed
cat FDRdrought | \sed s/^\Scaffold_//g | awk -v OFS="\t" '{ if ($13 <= 0.05) { print $1,$3,$3,$13} }' > CMH_383650.bed #383650

#I also subselect the sites that make through the Bonferonni threshold. 
cat FDRdrought | \sed s/^\Scaffold_//g | awk -v OFS="\t" '{ if ($14 <= 0.05) { print $1,$3,$3,$14} }' > CMH_BON.bed  #1419 loci

```

Then I use bedtools intersect to find intersections between beds. 

```
#!/bin/bash
#SBATCH --job-name=bedtools
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --time=5:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100GB

#conda environment
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner

#load files
cd /scratch/midway2/rozennpineau/drought/compare_sites_commongarden_drought/drought/
gwasbed=gwas_all.bed
cmhbed=CMH_all.bed
bedtools intersect -a $gwasbed -b $cmhbed -wa -wb -f 0.99 -r > intersect_all_drought_cmh_gwas.bed
#-r 1 -f 1 requires that there is at least 1 bp match, with 100% of the regions matching

gwasbed=gwas_893.bed
cmhbed=CMH_383650.bed
bedtools intersect -a $gwasbed -b $cmhbed -wa -wb -f 0.99 -r > intersect_significant_drought_cmh_gwas.bed

#for Bonferronni loci
gwasbed=gwas_893.bed
cmhbed=CMH_BON.bed 
bedtools intersect -a $gwasbed -b $cmhbed > intersect_BON_drought_cmh_gwas.bed #no intersection found

```


I further look for where in the genome the hits are : are there in genes ? I look into the gff file for this.
```
#prep bed file
awk '{OFS="\t";print $1, $2, $2}' common_significant_cmh_015.bed > cmh_015.bed
tail -n 13602 cmh_015.bed > cmh_015_noheader.bed

#add “Scaffold_” to overlap bed
awk -v OFS="\t" '{$1 = "Scaffold_" $1}1' cmh_015_noheader.bed > cmh_015_noheader_scaffold.bed
awk -v OFS="\t" '{$1 = "Scaffold_" $1}1' cmh_gwas_35.bed > cmh_gwas_35_scaffold.bed
awk -v OFS="\t" '{print $1, $2, $2}' common_significant_cmh_gwas_scaffold.bed > common_significant_cmh_gwas_scaffold_clean.bed
```

#bedtools intersect with gff
```
#remove header from bed file for bedtools intersect
bedtools intersect -a cmh_015_noheader_scaffold.bed -b /project/kreiner/data/genome/Atub_193_hap2.all.sorted.gff -wo > overlap_cmh_015_gff.bed
bedtools intersect -a cmh_gwas_35_scaffold_noheader.bed -b /project/kreiner/data/genome/Atub_193_hap2.all.sorted.gff -wo > overlap_gwas_cmh_35_gff.bed
bedtools intersect -a common_significant_cmh_gwas_scaffold_clean2.bed -b /project/kreiner/data/genome/Atub_193_hap2.all.sorted.gff -wo > overlap_gwas_cmh_gff.bed
```

#extract gene names
```
cut -f 12 overlap_cmh_015_gff.bed | grep "Similar" | cut -d" " -f3 | cut -d":" -f1 | sort | uniq > gwas_cmh_genes_35.txt
cut -f 12 overlap_gwas_cmh_35_gff.bed | grep "Similar" | cut -d" " -f3 | cut -d":" -f1 | sort | uniq > gwas_cmh_genes_35.txt
cut -f 12 overlap_gwas_cmh_gff.bed | grep "Similar" | cut -d" " -f3 | cut -d":" -f1 | uniq > gwas_cmh_genes.txt
```

### QQplot and Manhattan plos for CMH scan.

```
cd /scratch/midway2/rozennpineau/drought/compare_sites_commongarden_drought/drought

#file is very big, randomly subset 200,000 lines
awk 'NR > 1' FDRdrought | shuf -n 200000 > FDRdrought_rdn_subset #remove header

#add header back on
head -n 1 FDRdrought >  header
cat header FDRdrought_rdn_subse > FDRdrought_rdn_subset
#200001 ;ines with header
# plot qqplot and Manhattan plots from the CMH scan file that has all the info : FDRdrought
```

/project/kreiner/pairedenv_commongarden/normalized_SNPsonly_vcfs/drought_commongarden_merged_normalized_filtsnps_nomiss30.vcf.gz
#vcf file with the drought individuals 
/cds3/kreiner/drought/vcf_dont_touch/merged_numericChr.vcf.gz

### Thin result of CMH scan to only keep variants that are not in LD
Out of the 13,000+ variants, how many are actually independent ?
Steps: (1) filter vcf based on bed file, (2) make plink family files, (3) thin

(1) filter vcf based on bed file
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
vcftools --gzvcf drought_merged_numericChr.vcf.gz --bed cmh_015_scaffold.bed --out cmh_015 --recode

#note : check if you need to use the character (Scaffold_) and not numeric names!
#note : need for a slurm job, the terminal without a batch job was very slow!
```

Took ~30min. 
I have more than 13,603 variants, maybe several variant options at the same position ?


(2) make plink files
```
module load plink
vcf=cmh_015.recode.vcf

plink --vcf $vcf --out cmh_015_tothin --allow-extra-chr --recode --double-id 

```


This was going to work but the ID field in the vcf is missing, adding it now, and making it match with the association file, I need (1) Scaffold_ and (2) SNP changed to ID in the $assoc file

```
module load htslib
bgzip cmh_015.recode.vcf
bcftools tabix cmh_015.recode.vcf.gz
bcftools annotate --set-id +'%CHROM:%POS' cmh_015.recode.vcf.gz -o cmh_015_ID.vcf.gz

#plink family files
vcf=cmh_015_ID.vcf.gz
plink --vcf $vcf --out cmh_015_tothin --allow-extra-chr --recode --double-id 
```

(3) Clump 

```
#remove Scaffold from chromosome name in association file
sed 's/Scaffold_//g' $assoc > FDRdrought_numeric
#replace SNP with ID in association file
sed 's/SNP/ID/g' FDRdrought_numeric > FDRdrought_numeric_ID

assoc=FDRdrought_numeric_ID
plink --file cmh_015_tothin --clump $assoc --clump-p1 0.015 --clump-field FDR_p --clump-kb 100 --out cmh_015_clumped_100kb --allow-no-sex --allow-extra-chr --clump-snp-field ID

#output
Options in effect:
  --allow-extra-chr
  --allow-no-sex
  --clump FDRdrought_numeric_ID
  --clump-field FDR_p
  --clump-kb 100
  --clump-p1 0.015
  --clump-snp-field ID
  --file cmh_015_tothin
  --out cmh_015_clumped_100kb

257091 MB RAM detected; reserving 128545 MB for main workspace.
Possibly irregular .ped line.  Restarting scan, assuming multichar alleles.
.ped scan complete (for binary autoconversion).
Performing single-pass .bed write (14096 variants, 282 people).
--file: cmh_015_clumped_100kb-temporary.bed +
cmh_015_clumped_100kb-temporary.bim + cmh_015_clumped_100kb-temporary.fam
written.
14096 variants loaded from .bim file.
282 people (0 males, 0 females, 282 ambiguous) loaded from .fam.
Ambiguous sex IDs written to cmh_015_clumped_100kb.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 282 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.986247.
14096 variants and 282 people pass filters and QC.
Note: No phenotypes present.
Warning: '6:16669558' is missing from the main dataset, and is a top variant.
Warning: '8:10035631' is missing from the main dataset, and is a top variant.
Warning: '9:17730933' is missing from the main dataset, and is a top variant.
116748 more top variant IDs missing; see log file.
--clump: 4248 clumps formed from 13652 top variants.
Results written to cmh_015_clumped_100kb.clumped
```
**4248 clumps formed from 13652 top variants**


### How many of these clumps overlap with the drought loci?
```
#Make bed file from CMH clumps
awk -v OFS="\t" '{print $1, $4, $4}' cmh_015_clumped_100kb.clumped > cmh_015_clumped_100kb.clumped.bed
head -n 4249 cmh_015_clumped_100kb.clumped.bed > cmh_015_clumped_100kb_cleanclumped.bed #remove last lines that are just tabs
tail -n 4248 cmh_015_clumped_100kb_cleanclumped.bed > noheader #remove current header
cat header noheader > cmh_015_clumped_100kb_cleanclumped.bed #adjust header for bcftools to work
#the headers are not working, remove them ?

#Intersect with drought
cd /scratch/midway2/rozennpineau/drought/compare_sites_commongarden_drought/drought
gwasbed=/scratch/midway2/rozennpineau/drought/two_pulse_flexible_prop_2/drought_adapted_43clumps_noheader.bed
cmhbed=cmh_015_clumped_100kb_cleanclumped_noheader.bed
bedtools intersect -a $gwasbed -b $cmhbed -wa -wb -f 0.99 -r > intersect_clumped_drought_cmh_gwas.bed
bedtools intersect -a $gwasbed -b $cmhbed > intersect_clumped_drought_cmh_gwas.bed
#-r 1 -f 1 requires that there is at least 1 bp match, with 100% of the regions matching


```
