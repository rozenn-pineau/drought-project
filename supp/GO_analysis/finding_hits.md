### Identifying the genes with the lowest p-values from our GWAS

Filtering the association file to keep the p-values associated with each position. 
```
awk '$15 < 0.05' ancestry_corrected_inflated_gemma_gwas.assoc.txt > ancestry_corrected_inflated_gemma_gwas_893.assoc.txt

#keep rs, FDR and inflated p values
awk -F'\t' -vOFS='\t' '{print $13, $14, $14, $15, $16, $17, $18}' ancestry_corrected_inflated_gemma_gwas_893.assoc.txt > tmp

#Change chromosome to Scaffold
awk -F'\t' -vOFS='\t' '{ $1 = "Scaffold_" $1}1' tmp > tmp2

#header is "chrom   pos   pos2  FDR     bon     qlogp   inflation_p"
cat header_assoc tmp2 > ancestry_corrected_inflated_gemma_gwas_893_header.assoc.bed
```

Intersect with gff file and keep associated p-value.

```
bed=ancestry_corrected_inflated_gemma_gwas_893_header.assoc.bed
gff=/project/kreiner/data/genome/Atub_193_hap2.all.sorted.gff
bedtools intersect -b $bed -a $gff -wa -wb > intersect_893_gff.txt

```
Check p-value of lead snps, sorted based on inflation pvalue:

```
grep DDB1 intersect_893_gff.txt | sort -t$'\t' -k16,16 -g | head
#-g is for handling scientific notation
grep PER17 intersect_893_gff.txt | sort -t$'\t' -k16,16 -g | head
grep RKD5 intersect_893_gff.txt | sort -t$'\t' -k16,16 -g | head
grep TDC2 intersect_893_gff.txt | sort -t$'\t' -k16,16 -g | head
grep APSR1 intersect_893_gff.txt | sort -t$'\t' -k16,16 -g | head
```
