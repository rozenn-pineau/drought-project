### Randomization analysis for CMH scan / drought variants taking into account LD

** How many CMH scan SNPs overlap with the ancestry tracks in which the 43 clumps lead SNP for each are? 

(1) extract bed file that has start and end pos for each 43 clumps:

```
#add snpID info
awk 'BEGIN{OFS="\t"} NR==1{print $0,"snpID"; next} {print $0, NR-1}' drought_adapted_43clumps.bed > drought_adapted_43clumps_snpID.bed


anc_len=/scratch/midway3/rozennpineau/drought/randomization/ancestry_block_len_genome_tab.txt
bed=/scratch/midway3/rozennpineau/drought/randomization/1_overlap_cmh_gwas_observed/drought_adapted_43clumps_snpID.bed

bedtools intersect -a $anc_len -b $bed -wa -wb \
| awk 'BEGIN{OFS="\t"} {print $1,$2,$3, $4,$5,$(NF)}' > intersect_43drought_blocklen_tab.txt
#entries from the ancestry block file that I want to keep that intersect the 43 drought adapted loci.
#information about the snpID in the last column

#add header
cat header intersect_43drought_blocklen_tab.txt > intersect_43drought_blocklen_tab_header.txt

```
**Note : do we want to keep SNPs that are at the edge between two ancestry tracts ?
--> maybe I need to re-work the script that calculates the lengths of each tract to make sure it starts at index + 1 for the next tract

(2) count the number of CMH SNPs that fall within those ancestry tracts. This means I have a number per SNP per sample.

```

block43=intersect_43drought_blocklen_tab_header.txt
cmhbed=/scratch/midway2/rozennpineau/drought/compare_sites_commongarden_drought/drought/CMH_FDR01.bed

bedtools intersect -a $block43 -b $cmhbed -c > intersect_cmh_43block.txt
#-c For each feature in A, report the number of overlapping features in B

#then look for length of overlaps
gwasbed=drought_adapted_43clumps.bed


```
#tab delimitation for /scratch/midway3/rozennpineau/drought/randomization/ancestry_block_len_genome.txt
awk '{print $1, $2, $3, $4, $5}' OFS='\t' ancestry_block_len_genome.txt > ancestry_block_len_genome_tab.txt


[calculate_ancestry_block_lengths_v2.R](https://github.com/rozenn-pineau/drought-project/tree/main/cmh_gwas_randomization/calculate_ancestry_block_lengths_v2.R)
