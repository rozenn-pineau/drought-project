### andomization analysis for CMH scan / drought variants taking into account LD

** How many CMH scan SNPs overlap with the ancestry tracks in which the 43 clumps lead SNP for each are? 

(1) extract bed file that has start and end pos for each 43 clumps : [calculate_ancestry_block_lengths_v2.R](https://github.com/rozenn-pineau/drought-project/tree/main/cmh_gwas_randomization/calculate_ancestry_block_lengths_v2.R)

(2) count the number of CMH SNPs that fall within those ancestry tracts.

#then look for length of overlaps
bed=/scratch/midway2/rozennpineau/drought/two_pulse_flexible_prop_2/clump_summary_scaffold_names.bed
gff=/scratch/midway2/rozennpineau/drought/two_pulse_flexible_prop_2/intersect_43clumps_gff_gene.txt
bedtools intersect -a $bed -b $gff -c
#-c For each feature in A, report the number of overlapping features in B
