## Admixture mapping

### Using GEMMA with the relatedness matrix

### Prepping the data for GWAS

To use the ancestry calls with gemma, I first had to reformat it to bimbam format :

```

#!/bin/bash


cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2/NAs
file=/scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2/NAs/two_pulse_flexible_prop_2_values_cutoff_08_v2.txt
output=/scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2/NAs/two_pulse_flexible_prop_2_values_cutoff_08_v2.gemma
#remove first line
head -n +1 $file > samp.order

tail -n -786261 $file > tmp

#remove outlier columns
awk -F'\t' '{for(i=1;i<=NF;i++) {if($i == "P16_Nat_1_T") printf("Column %d is outlier #1\n", i-1)}}' samp.order
awk -F'\t' '{for(i=1;i<=NF;i++) {if($i == "P12_Nat_14_T") printf("Column %d is outlier #2\n", i-1)}}' samp.order
#columns 29 and 104
cut -d$'\t' --complement -f 29,104 tmp > tmp1

#add fake ref to alt allele columns in positions 3 and 4
awk 'BEGIN { FS = OFS = "\t" } {
    for (i = 1; i <= NF; i++) {
        if (i == 3) printf "A\tT\t"; # Add "A" and "T" before column 3
        printf "%s%s", $i, (i < NF ? OFS : ""); # Print original column with tab separator
    }
    print ""; # Add newline after each row
}' tmp1 > tmp2

#replace first tab by ":", replace every remaining tab by commas, and delete column 2
awk -F'\t' '{$1 = $1 ":" $2; print}' tmp2 > tmp3
awk 'BEGIN { FS=" "; OFS="," } {$1=$1; print}' tmp3 > tmp4
cut -d ',' --complement -f2 tmp4 > $output


rm -I tmp*

```

Running GEMMA (we allowed up to 20% missing data per sample for each SNP):

```
#!/bin/bash
#SBATCH --job-name=gemma
#SBATCH --output=gemma.out
#SBATCH --error=gemma.err
#SBATCH --time=25:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100G   # memory per cpu-core

cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/gemma_gwas/two_pulse_flexible_prop_2

anc_file=/scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2/NAs/two_pulse_flexible_prop_2_values_cutoff_08_v2.gemma
phenos=/scratch/midway3/rozennpineau/drought/ancestry_hmm/manual_gwas/phenos_ordered_noout.txt
geno_relatedness=/scratch/midway3/rozennpineau/drought/ancestry_hmm/gemma_gwas/relatedness_matrix_geno_ordered.txt

module load gemma
gemma -g $anc_file -p $phenos -miss 0.1 -lm #gwas not corrected for population structure
#-miss 0.1 is to accept up to 10% missing data per sample per SNP, will be interpolated from existing data

#gemma gwas corrected for pop structure
##with genotype-based relatedness matrix
gemma -g $anc_file -p $phenos -k $geno_relatedness -km 1 -miss 0.1 -lmm -o geno_corrected_gemma_gwas

##with anc covariate, first generating relatedness matrix (-gk)

gemma -g $anc_file -p $phenos -gk 1 -o ancestryrelatedness #generate relatedness
gemma -g $anc_file  -p $phenos -k output/ancestryrelatedness.cXX.txt -miss 0.1 -lmm -o ancestry_corrected_gemma_gwas

```

For the genotype-based relatedness matrix, a bit of sample reordering was necessary and done in R :

```
#goal: to re order the relatedness matrix built from genotypes such that it matches the ancestry call order
rm(list= ls())

#load relatedness matrix 
mat <- read.table("/project/kreiner/drought/gwas_episode3/output/merged_numericChr_nooutliers_maf01.cXX.txt", sep = "\t", header = F)#280 x 280

#load order of samples in rm
order_mat <- read.table("/project/kreiner/drought/gwas_episode3/prep_files/merged_numericChr_nooutliers_maf01.fam", sep = " ", header=F)

#load order of samples in ancestry calls
order_goal <- read.table("/scratch/midway3/rozennpineau/drought/ancestry_hmm/manual_gwas/samp_order_noout.txt", sep = "\t", header = F)

#order matrix based on order of samples in order_goal
#template : vector1[order(match(vector1,vector2))] 
new_order <- order(match(order_mat[,1], order_goal[,1]))
new_mat <- mat[new_order, new_order]

#export
write.table(new_mat, "/scratch/midway3/rozennpineau/drought/ancestry_hmm/gemma_gwas/relatedness_matrix_geno_ordered.txt", sep = "\t", col.names = F, row.names = F)
```
## Analyzing associated sites in R and correcting p-values

We analyzed the runs from different probability cutoffs (0.5 to 0.9) and different NA content tolerance in GEMMA command line (0.05 -default-, 0.1 and 0.2). 
We compared QQplots and Manhanttan plots from GWAS with and without covariate (ancestry and genotype-based relatedness matrices). 
Plot for Figure 3 admixture mapping analyses : [figure3.R](https://github.com/rozenn-pineau/drought-project/blob/main/figure3.R)

## Clumping based on LD
Next step is to clump together the sites that may be in LD. We do this using plink.

First, we export the **inflated** and **FDR corrected** p-values from the association file (see above script). We will use these inflated p values to clump sites based on LD (and significance).




(1) make vcf based on table
```
#!/bin/bash
#SBATCH --job-name=vcf_plink
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --time=20:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=50GB

cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2/

# Input and output files
input_file="two_pulse_flexible_prop_2_values.txt"
output_file="two_pulse_flexible_prop_2_values.vcf"

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
awk 'NR <= 5 {print; next} {OFS="\t"; $3 = $1 ":" $2; print}' two_pulse_flexible_prop_2_values.vcf > two_pulse_flexible_prop_2_values_ID.vcf
```

(2) create plink family files (.ped and .map)
```
cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2
#make plink family files and create ID based on position information
plink --vcf two_pulse_flexible_prop_2_values_ID.vcf --out two_pulse_flexible_prop_2_values --allow-extra-chr --recode --double-id 
```

(3) clump sites using plink

```
#add ID column to association file
sed -e 's/rs/ID/g' ancestry_corrected_inflated_gemma_gwas.assoc.txt > ancestry_corrected_inflated_gemma_gwas_ID.assoc.txt        

#run plink clump
plink --file two_pulse_flexible_prop_2_values --clump ancestry_corrected_inflated_gemma_gwas_ID.assoc.txt --clump-p1 6.5830432740603e-05 --clump-field p_wald --clump-kb 100 --out two_pulse_flexible_prop_2_clumped_100kb_43_clumps --allow-no-sex --allow-extra-chr --clump-snp-field ID

#value for pvalue chosen by selecting the closest value to 0.05 from the transformed (FDR) vector
#bfile is 0.9 cutoff for ancestry calls
#association file from gemma gwas with --miss 0.2
```

Output from plink --clump :


```
Options in effect:
  --allow-extra-chr
  --allow-no-sex
  --clump ancestry_corrected_inflated_gemma_gwas_ID.assoc.txt
  --clump-field p_wald
  --clump-kb 100
  --clump-p1 6.5830432740603e-05
  --clump-snp-field ID
  --file two_pulse_flexible_prop_2_values
  --out two_pulse_flexible_prop_2_clumped_100kb_43_clumps

257091 MB RAM detected; reserving 128545 MB for main workspace.
.ped scan complete (for binary autoconversion).
Performing single-pass .bed write (786261 variants, 282 people).
--file: two_pulse_flexible_prop_2_clumped_100kb_43_clumps-temporary.bed +
two_pulse_flexible_prop_2_clumped_100kb_43_clumps-temporary.bim +
two_pulse_flexible_prop_2_clumped_100kb_43_clumps-temporary.fam written.
786261 variants loaded from .bim file.
282 people (0 males, 0 females, 282 ambiguous) loaded from .fam.
Ambiguous sex IDs written to
two_pulse_flexible_prop_2_clumped_100kb_43_clumps.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 282 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.889586.
786261 variants and 282 people pass filters and QC.
Note: No phenotypes present.
--clump: 43 clumps formed from 1032 top variants.
Results written to two_pulse_flexible_prop_2_clumped_100kb_43_clumps.clumped .

```
43 clumps !!


Prep the files to analyze allelic trajectories :

(1) make bed file to filter ancestry calls vcf file

```
#!/bin/bash

file=two_pulse_flexible_prop_2_clumped_100kb_43_clumps.clumped

awk -v OFS='\t' '{ print $1, $4, $4, $5}' $file > drought_adapted_43clumps.bed #43 loci in this file

```

(2) filter the vcf file based on the bed file

```
module load vcftools
cd /scratch/midway3/rozennpineau/drought/ancestry_hmm/run_full_genome/two_pulse_flexible_prop_2
vcftools --vcf two_pulse_flexible_prop_2_values_ID.vcf --bed drought_adapted_43clumps.bed --out drought_adapted_43clumps --recode


#output 
Parameters as interpreted:
	--vcf two_pulse_flexible_prop_2_values_ID.vcf
	--out drought_adapted_43clumps
	--recode
	--bed drought_adapted_43clumps.bed

After filtering, kept 282 out of 282 Individuals
Outputting VCF file...
	Read 44 BED file entries.
After filtering, kept 43 out of a possible 786261 Sites
Run Time = 3.00 seconds
````
43 loci in the output vcf file.


Getting the ancestry calls (in the form of genotypes) from the filtered vcf for more downstream analyses :

```
module load htslib
bgzip -f drought_adapted_43clumps.recode.vcf
tabix -f drought_adapted_43clumps.recode.vcf.gz
act-conda
bcftools query -f '%CHROM %POS  %REF  %ALT [ %GT]\n' drought_adapted_43clumps.recode.vcf.gz > drought_adapted_43clumps_GT.txt #43 loci

#download on laptop for further analyses in R.
scp rozennpineau@midway3.rcc.uchicago.edu:/scratch/midway2/rozennpineau/drought/two_pulse_flexible_prop_2/drought_adapted_43clumps_GT.txt /Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My\ Drive/Work/9.Science/1.DroughtProject/1.analyses/data/6.trajectories/
```


### Checking : Does ancestry predict response to drought ?

Our expectation is that var. rudis ancestry is better adapted to drought than var. tuberculatus. 
Do we see this in our ancestry calls ? This is also simply a way to make sure that our pipeline was coded correctly. 

