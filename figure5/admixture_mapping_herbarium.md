## Admixture mapping on the herbarium dataset. 

### Step (1) : keep the variants common to the ancestry variant file and the herbarium variant file
To build the file that will be used to call ancestry at each of the herbarium site, we keep the intersection between the ancestry variant file and the herbairum variant file using bcftools isec.

```
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner

#bgzip and tabix the vcf files

anc=/scratch/midway3/rozennpineau/drought/ancestry_hmm/prep_ancestry/fst/ancestry_fst_q75.vcf.gz

herb=/scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/merged_numericChr.vcf.gz

cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_more_sites/

bcftools isec -p herb_ancestry $anc $herb

```

At this step, we have two variant files with the same variants for the two populations (976662 sites). I extract the position and chromosome information (bed) to calculate rho between each site. 

### Step (2) : calculate rho between each site

Extract chromosome and position information from the vcfs :

```
awk 'BEGIN {OFS="\t"} !/^#/ {print $1, $2-1, $2}' ancestry_common_976662.vcf  > ancestry_common_976662.be

```

Rscript to calculate rho between sites : [calculate_ldhat_between_sites.Rmd](https://github.com/rozenn-pineau/Drought-paper/blob/main/calculate_ldhat_between_sites.Rmd).

Because we exclude the sites that are outside of the region defined by the recombination map (the monotonic spline does not extrapolate outside of boundaries), we kept track of which sites to filter out of the variant files to further filter the vcf files.

```
#make bed file
awk '{print $1,$2,$3}' ancestry_herb_common_974134.ld | tail -n -974134 > ancestry_herb_common_974134.bed

#filter
module load vcftools
bed=/scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_more_sites/ancestry_herb_common_974134.bed
vcf1=/scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_more_sites/herb_ancestry/herb_common_976662.vcf
vcftools --vcf $vcf1 --bed $bed --out herb_common_974134.vcf --recode

vcf2=/scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_more_sites/herb_ancestry/ancestry_common_976662.vcf
vcftools --vcf $vcf2 --bed $bed --out ancestry_common_974134.vcf --recode

```

### Step (3) : split ancestry vcf into tuberculatus versus rudis files

The ancestry file both has var rudis and var tuberculatus samples, that we will now split into two :

```
bcftools view -S var_rudis_samp.txt ancestry_common_974133.vcf > rudis_common_974133.vcf
bcftools view -S var_tub_samp.txt ancestry_common_974133.vcf > tub_common_974133.vcf

#check # of variants and samples
#44 samples for rudis and 21 samples for tuberculatus
#974133 variants for both
```


### Step (4) : extract allele counts from the ancestry variant files

[genotype_to_allele_counts.awk](https://github.com/rozenn-pineau/Drought-paper/blob/main/genotype_to_allele_counts.awk)

```
/scratch/midway3/rozennpineau/drought/scripts/genotype_to_allele_counts.awk tub_common_974133.vcf > tub_common_974133.allele_counts

/scratch/midway3/rozennpineau/drought/scripts/genotype_to_allele_counts.awk rudis_common_974133.vcf > rudis_common_974133.allele_counts

```

### Step (5) : get genotypes for drought panel

2,0 for homs reference
1,1 for hets
0,2 for homs alternative

[vcf_to_read_counts.awk](https://github.com/rozenn-pineau/Drought-paper/blob/main/vcf_to_read_counts.awk)

### Step (6) : put the file together

Paste columns together to make the full file, following instructions on ancestry_hmm github page. 

```


#var rudis allele counts 
awk '{OFS="\t"; print $3,$4}' tub_common_974133.allele_counts > var_rud.allele_counts

#rho column
awk '{OFS="\t"; print $4}' ancestry_herb_common_974133.ld > ancestry_herb_common_974133_rho.ld

#assemble chrom, pos, var tub allele counts, var rudis allele counts, rho, then sample read counts
paste tub_common_974133.allele_counts var_rud.allele_counts ancestry_herb_common_974133_rho.ld herb_common_974133.read_counts > herb_common_974133_input_file.txt


```

### Step (7) : make sample file


### Process output from ancestry_hmm

```
folder_name=$(basename "$PWD")
output_file="${folder_name}_values_09.txt"
cutoff=0.9

#initialize dataset with chrom and pos
awk 'BEGIN { OFS="\t" } {print $1,$2}' HB0900.posterior > $output_file

#loop through each individual and paste information to previous version of the file
for file in *.posterior; do

    header_name=$(basename "$file" .posterior)

    awk -v header="$header_name" 'NR == 1 { print header; }
    
    NR > 1 {
    
    result = "NA";
    if ($3 > $cutoff) {result = "0"}
    else if ($4 > $cutoff) {result = "1"}
    else if ($5 > $cutoff) {result = "1"}
    else if ($6 > $cutoff) {result = "2"}
    else if ($7 > $cutoff) {result = "2"}
    else if ($8 > $cutoff) {result = "2"}
    print result;
    }' $file > tmp

    paste $output_file tmp > tmp2

    mv tmp2 $output_file

done
```

This is for every site in common between the herbarium variant file and the ancestry panels. We are interest in the trajectories of drought-informative alleles. Let's filter the file for those sites only.  

### (1) filter ancestry call table for drought informative sites

```
cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_more_sites/1_two_pulse_flexible_proportions

# Input files
bed=/scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_more_sites/FDR_non_clumped_significant_sites.bed
assoc=/scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_more_sites/1_two_pulse_flexible_proportions/anc_calls_herbarium_all_sites_09.txt
out=/scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_more_sites/1_two_pulse_flexible_proportions/anc_calls_herbarium_893sites_09.txt

# Define column positions for chrom and pos in the TXT file (1-based index)
CHROM_COL=1
POS_COL=2

# Convert column positions to AWK's 1-based indexing
awk -v chrom_col="$CHROM_COL" -v pos_col="$POS_COL" '
    NR==FNR {sites[$1][$3]; next}
    {
        chrom = $chrom_col;
        pos = $pos_col;
        if (chrom in sites)
            for (s in sites[chrom])
                if (pos >= s && pos <= s)
                    print $0;
    }' "$bed" "$assoc" > "$out"

echo "Filtering complete. Results saved in $out"

mv anc_calls_herbarium_893sites_09.txt anc_calls_herbarium_680sites_09.txt
#680 sites in common

#add header back to ancestry call file
```
We now have 680 sites with ancestry calls, that are drought informative sites. 

### (2) Filter association file for these 680 drought informative sites
```
#make bed from vcf file
vcf=anc_calls_herbarium_680sites_09.vcf 
bed=anc_calls_herbarium_680sites_09.bed

awk 'BEGIN {OFS="\t"} 
     !/^#/ {print $1, $2-1, $2}' "$vcf" > "$bed" #position in vcf is shifted by one

# Input files
bed=anc_calls_herbarium_680sites_09.bed
assoc=/scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2/ancestry_corrected_inflated_gemma_gwas_ID.assoc.txt
out=FDR_non_clumped_680sites.assoc.txt

# Define column positions for chrom and pos in the TXT file (1-based index)
CHROM_COL=13  
POS_COL=14    


# Convert column positions to AWK's 1-based indexing
awk -v chrom_col="$CHROM_COL" -v pos_col="$POS_COL" '
    NR==FNR {sites[$1][$3]; next}  
    { 
        chrom = $chrom_col; 
        pos = $pos_col;
        if (chrom in sites) 
            for (s in sites[chrom]) 
                if (pos >= s && pos <= s) 
                    print $0;
    }' "$bed" "$assoc" > "$out"

echo "Filtering complete. Results saved in $out"

```
### (3) make vcf file for plink LD thinning and drop variants with more than 75% missing data

**make vcf file**

Plink takes in a vcf file, that we create based on the ancestry call table. 

```
cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/herbarium/2_more_sites/1_two_pulse_flexible_proportions/

# Input and output files
input_file="anc_calls_herbarium_680sites_09.txt"
output_file="anc_calls_herbarium_680sites_09.vcf"

# Create the VCF header
cat <<EOL > $output_file
##fileformat=VCFv4.2
##source=CustomScript
##INFO=<ID=.,Number=1,Type=String,Description="Custom info">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
EOL
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$(head -1 $input_file | cut -f3- | tr -s ' ' '\t')" >> $output_file


# Transform the input data into VCF format
tail -n +2 $input_file | while read -r line; do
    # Parse the fields
    chrom=$(echo "$line" | awk '{print $1}')
    pos=$(echo "$line" | awk '{print $2}')
    genotypes=$(echo "$line" | cut -f3-)

    # Convert genotypes to VCF GT format
    formatted_genotypes=$(echo "$genotypes" | awk '
    {
        for (i = 1; i <= NF; i++) {
            if ($i == 0) $i = "0|0";      # Homozygous reference
            else if ($i == 1) $i = "1|0";      # Heterozygote
            else if ($i == 2) $i = "1|1"; # Homozygous alternate
            else $i = "./.";              # Missing data or unrecognized value
        }
        print $0
    }' | sed 's/ /\t/g')

    # Placeholder values for REF, ALT, ID, QUAL, FILTER, INFO
    ref="A" # Adjust this to your data's REF allele
    alt="T" # Adjust this to your data's ALT allele
    id="." # No variant ID provided
    qual="."
    filter="PASS"
    info="."
    format="GT"

    # Print the VCF record
    echo -e "$chrom\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$format\t$formatted_genotypes" >> $output_file
done

echo "VCF file created: $output_file"


#fill in the ID field by combining chrom and pos
awk 'NR <= 5 {print; next} {OFS="\t"; $3 = $1 ":" $2; print}' anc_calls_herbarium_680sites_09.vcf > anc_calls_herbarium_680sites_09_ID.vcf
```

**drop loci with more than 75% missing data**


### (4) Clumping based on LD

```
assoc=FDR_non_clumped_680sites_header.assoc.txt

#without missing loci threshold

plink --file anc_calls_herbarium_680sites_09_ID --clump $assoc --clump-p1 6.5830432740603e-05 --clump-field p_wald --clump-kb 100 --out two_pulse_flexible_prop_2_clumped_100kb --allow-no-sex --allow-extra-chr --clump-snp-field ID

#output
Logging to two_pulse_flexible_prop_2_clumped_100kb.log.
Options in effect:
  --allow-extra-chr
  --allow-no-sex
  --clump FDR_non_clumped_680sites_header.assoc.txt
  --clump-field p_wald
  --clump-kb 100
  --clump-p1 6.5830432740603e-05
  --clump-snp-field ID
  --file anc_calls_herbarium_680sites_09_ID
  --out two_pulse_flexible_prop_2_clumped_100kb

257091 MB RAM detected; reserving 128545 MB for main workspace.
.ped scan complete (for binary autoconversion).
Performing single-pass .bed write (680 variants, 108 people).
--file: two_pulse_flexible_prop_2_clumped_100kb-temporary.bed +
two_pulse_flexible_prop_2_clumped_100kb-temporary.bim +
two_pulse_flexible_prop_2_clumped_100kb-temporary.fam written.
680 variants loaded from .bim file.
108 people (0 males, 0 females, 108 ambiguous) loaded from .fam.
Ambiguous sex IDs written to two_pulse_flexible_prop_2_clumped_100kb.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 108 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.677124.
680 variants and 108 people pass filters and QC.
Note: No phenotypes present.
--clump: 33 clumps formed from 680 top variants.
Results written to two_pulse_flexible_prop_2_clumped_100kb.clumped .

#--geno <threshold>: removes SNPs where the missing genotype proportion exceeds the provided threshold. --geno 0.25 removes loci with more than 25% missing genotype data, which corresponds to keeping loci with at least 75% non-missing genotype data.

```


### (5) extract clumped variants 
(1) make bed file to filter ancestry calls vcf file
```
awk -v OFS='\t' '{ print $1, $4, $4, $5}' two_pulse_flexible_prop_2_clumped_100kb.clumped > herbarium_filtered_sites.bed
```
(2) filter the vcf file based on the bed file
```
module load vcftools
vcftools --vcf anc_calls_herbarium_680sites_09_ID.vcf --bed herbarium_filtered_sites.bed --out anc_calls_herbarium_33sites --recode

Parameters as interpreted:
	--vcf anc_calls_herbarium_680sites_09_ID.vcf
	--out anc_calls_herbarium_33sites
	--recode
	--bed herbarium_filtered_sites.bed

After filtering, kept 108 out of 108 Individuals
Outputting VCF file...
	Read 34 BED file entries.
After filtering, kept 33 out of a possible 680 Sites
Run Time = 0.00 seconds
```
(3) getting ancestry calls for downstream analyses
```
module load htslib
bgzip -f anc_calls_herbarium_33sites.recode.vcf 
tabix -f anc_calls_herbarium_33sites.recode.vcf.gz

bcftools query -f '%CHROM %POS  %REF  %ALT [ %GT]\n' anc_calls_herbarium_33sites.recode.vcf.gz > anc_calls_herbarium_33sites_GT.txt
```

