### GO enrichment analyses and sites of potential interest

We identified 36 sites significantly associated with drought adaptation. To test whether those sites are in functional regions of the genome with known effects, we extracted the regions from the annotated genome file.

The gff file has "Scaffold_" as chromosome names, so I need to update the bed chromosome names :

```
awk -F'\t' -vOFS='\t' '{ $1 = "Scaffold_" $1}1' ancestry_gwas_filtered_sites.bed > ancestry_gwas_filtered_sites_scaffold_names.bed 
```
Intersect the bed file with the gff file using bedtools :

```
bedtools intersect -b ancestry_gwas_filtered_sites_scaffold_names.bed -a /project/kreiner/data/genome/Atub_193_hap2.all.sorted.gff > intersect_FDR_gff_enriched_genes.txt
```
Extract the GO terms and Note field from the file : 

```
cut -d ";" -f 9 intersect_FDR_gff_enriched_genes.txt | grep GO > intersect_FDR_gff_GO_terms_bon.txt
cut -d ";" -f 10 intersect_FDR_gff_enriched_genes.txt | grep Note > intersect_FDR_gff_note.txt
```


### snpEff 

I followed the below steps to find the effect of the mutation on the protein function (thank you, Jake!) : 

```
#make a conda environment for snpeff
conda create -n snpeff

#activate that environment
conda activate snpeff

#install snpeff
conda install bioconda::snpeff

#move to the directory where conda installed snpeff
cd .conda/envs/snpeff/share/snpeff-5.2-1/

#make a new directory to store your database(s)
mkdir data

#change into there and make a directory for your database
cd data
mkdir Atub_193_hap2

#change into there and copy the genome sequence to a file called sequences.fa
cd Atub_193_hap2
cp /project/kreiner/data/genome/Atub_193_hap2.fasta sequences.fa

#either install gffread in this environment, or in my case, I use a different conda environment
#conda deactivate
#conda create -n gffread
#conda activate gffread
conda install bioconda::gffread

#   .. loaded 35100 genomic features from /project/kreiner/data/genome/Atub_193_hap2.all.sorted.gff

#use gffread to convert the gff file into gtf
gffread -E /project/kreiner/data/genome/Atub_193_hap2.all.sorted.gff -T -o genes.gtf

#use gffread to create cds and protein sequence files
gffread -x cds.fa -g /project/kreiner/data/genome/Atub_193_hap2.fasta /project/kreiner/data/genome/Atub_193_hap2.all.sorted.gff
gffread -y protein.fa -g /project/kreiner/data/genome/Atub_193_hap2.fasta /project/kreiner/data/genome/Atub_193_hap2.all.sorted.gff

#move back over to your snpeff environment and add the database that you want to build to the snpeff config file
nano snpEff.config

# add the following lines to the config file right below the '# Non-standard Databases' lines, above 'Homo sapiens (hg19) (UCSC)'

# Atuberculatus genome, version 193_hap2
Atub_193_hap2.genome : Atub_193_hap2

#build the database
snpEff build -gtf22 -v Atub_193_hap2

#run snpeff on the vcf
#filter the vcf based on the bed file
bedtools intersect -header -b /scratch/midway2/rozennpineau/drought/compare_sites_commongarden_drought/FDR_significant_drought_sites.bed -a /scratch/midway2/rozennpineau/drought/two_pulse_flexible_prop_2/two_pulse_flexible_prop_2_values_ID.vcf > FDR_significant_drought_sites.vcf
#update scaffold names to include "Scaffold_"
awk -F'\t' -vOFS='\t' '{$1 = "Scaffold_" $1}1' FDR_significant_drought_sites.vcf > FDR_significant_drought_sites_scaffold_names.vcf
snpEff ann Atub_193_hap2 FDR_significant_drought_sites_scaffold_names.vcf > FDR_significant_drought_sites_ann.vcf

```


