# drought-project
Scripts and analyses for the drought paper.


# Outline
[Figure 1 analyses](#Figure-1-analyses)

[Trade-offs analyses](#Trade-offs-analyses)

- Herbicide resistance estimations


[Figure 2 analyses](#Figure-2-analyses)
  
- Estimating global ancestry using ADMIXTURE
  
[Figure 3 analyses](#Figure-3-analyses)

- Ancestry mapping using ancestry hmm on the contemporary dataset
  
[Figure 4 analyses](#Figure-4-analyses)

[Figure 5 analyses](#Figure-5-analyses)

[Figure 6 analyses](#Figure-6-analyses)


# Figure 1 analyses
Drought LD50 calculations, multiple linear regression models and Wilcoxon tests : [figure1_v2.R](https://github.com/rozenn-pineau/drought-project/blob/main/figure1_v2.R) and [Figure1E.Rmd](https://github.com/rozenn-pineau/drought-project/blob/main/Figure1E.Rmd).

# Trade-offs analyses
### Herbicide resistance estimations

To quantify resistance in drought samples, we calculated the depth at EPSPS position in each sample. 

We first extracted the lines from the gff files that have "EPSPS" to identify the gene location in the genome: 

```
#extract lines with EPSPS
grep "EPSPS" Atub_193_hap2.all.sorted.gff > /scratch/midway3/rozennpineau/drought/intersect_gff.txt

```
There was one site within the genome 16 scaffolds : 
```
Scaffold_6      maker   mRNA    17472580        17484501
```
We further used mosdepth to extract the depth at this site for every sample (from the bam files)  : 

```
cd /cds3/kreiner/drought/bams
bed=/scratch/midway3/rozennpineau/drought/EPSPS/EPSPS.bed
for bam in *.bam
do
mosdepth -n -x --by $bed /scratch/midway3/rozennpineau/drought/EPSPS/mosdepth_output/$bam $bam
done

#From the mosdepth output, we extracted the depth into one file with the corresponding sample name.
for file in *regions.bed.gz
do
zcat $file | awk 'BEGIN {OFS="\t"}; {print $4}' >> depth.txt
done

#add the file name info in column 1 (extract name from file name?)
ls *regions.bed.gz > file_names.tmp #make a file of file names
cut -d "_" -f 1,2 file_names.tmp > samp_ID.tmp #extract important IF info
paste samp_ID.tmp depth.txt > ID_depth.txt #tab delim file
```
After extracting te depths for each sample at EPSPS gene, we need to extract the genome wide mean coverage for normalization : 

```
#make bed file from gff
awk 'BEGIN {OFS="\t"}; {print $1,$4,$5}' /project/kreiner/data/genome/Atub_193_hap2.all.sorted.gff >gff.bed
uniq gff.bed > gff_clean.bed

#extract depth for each coding site
bed=/scratch/midway3/rozennpineau/drought/EPSPS/gff_clean.bed

cd /cds3/kreiner/drought/bams

for bam in *.bam
do
mosdepth -n -x --by $bed /scratch/midway3/rozennpineau/drought/EPSPS/gff_depths/$bam $bam 
done

#gather in one file  : paste them together 

zcat AT27_827_193_2.dd.bam.regions.bed.gz | awk 'BEGIN {OFS="\t"}; {print $1,$2,$3}' > depth.txt
#410700 lines

for file in *regions.bed.gz
do
zcat $file | awk 'BEGIN {OFS="\t"}; {print $4}' > tmp
cp depth.txt tmp2
paste tmp2 tmp > depth.txt
done

#same file order as the EPSPS depth file

```
Rscript to plot results from trade-offs analyses : [trade_offs_control.Rmd](https://github.com/rozenn-pineau/drought-project/blob/main/trade_offs_control.Rmd)

# Figure 3 analyses

### Estimating global ancestry using ADMIXTURE
To calculate genome-wide ancestry, we used ADMIXTURE version 1.3.0 (Alexander et al., 2009). 
We tested different values for K and settled on K=2: 

```
#!/bin/bash
#SBATCH --job-name=admixture
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --time=15:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100G   # memory per cpu-core

module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner

my_bed=/scratch/midway3/rozennpineau/drought/ancestry_hmm/prep_ancestry/commongarden_allfiltsnps_193_hap2_numericChr_filt.bed

for K in 2;

do admixture --cv $my_bed $K ;

done
```

# Figure 3 analyses

## Ancestry mapping using ancestry hmm on the contemporary dataset

[Ancestry_hmm](https://github.com/russcd/Ancestry_HMM) : tool to infer ancestry at input positions in the genome (Corbett-Detig, R. and Nielsen, R., 2017.)


### Step (1) : define var rudis versus and var tuberculatus pure ancestry individuals
To define the "ancestry panel" and not lose any of the drought data, we used an additional dataset, variants from common garden experiments used in Kreiner et al, 2022 that were realigned on a newer version of the reference genome. We ran structure with k=2, and kept individuals with pure ancestry with a threshold of 0.00001. We identified 44 pure var. rudis and 21 pure var. tuberculatus samples.

```
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner

my_bed=/scratch/midway3/rozennpineau/drought/ancestry_hmm/prep_ancestry/commongarden_allfiltsnps_193_hap2_numericChr_filt.bed

for K in 2;

do admixture --cv $my_bed $K ;

done
```

### Step (2) : calculate Fst on ancestry sites with MAF <0.05
We filtered the pure ancestry variant file (50,811,811 mutations) for minimum allele frequency of 0.05. To focus the analysis on ancestry-informative sites, we kept the sites with Fst values situated in the 75th percentile and above (2,089,620 variants, Fst calculated using VCFtools version 0.1.16).


```
#calculate fst per site with magf > 0.05
module load vcftools

VCF=/scratch/midway3/rozennpineau/drought/ancestry_hmm/prep_ancestry/ancestry.vcf.recode.vcf
vcftools --vcf ${VCF} --maf 0.05 --recode --stdout > fst/ancestry_maf.vcf

vcftools --vcf fst/ancestry_maf.vcf \
--weir-fst-pop var_rudis_samp.list \
--weir-fst-pop var_tub_samp.list \
--out fst/ancestry_maf

```

Rscript to choose Fst threshold : [fst_on_ancestry.Rmd](https://github.com/rozenn-pineau/drought-project/blob/main/fst_on_ancestry.Rmd)

### Step (3) : keep the variants common to the ancestry variant file and the drought variant file
The ancestry was filtered for the high Fst sites. This filtered file was used to keep the intersection between the ancestry variant file and the drought variant file (bcftools isec).
At this step, we have two variant files with the same variants for the two populations.
Now, we need to generate the rho information between each site.


### Step (4) : calculate rho between each site
Based on an estimation of ld for (nearly) the whole genome, we calculated, chromosome per chromosome, the linkage value between two consecutive sites. To do this, we fit a monotonic spline on the cumulative distribution of the mean recombination rates. 

Important note: we had to exclude the sites that were outside of the region defined by the recombination map (the monotonic spline does not extrapolate outside of boundaries). Then, we kept track of which sites to filter out of the variant files. 

Rscript to calculate rho between sites : [calculate_ldhat_between_sites.Rmd](https://github.com/rozenn-pineau/Drought-paper/blob/main/calculate_ldhat_between_sites.Rmd).


### Step (5) : extract allele counts from the ancestry variant file

[genotype_to_allele_counts.awk](https://github.com/rozenn-pineau/drought-project/blob/main/genotype_to_allele_counts.awk)

exanple output :
22 22 for 22 individuals, all heterozygotes
00 44 for 22 individuals, all homs alternative

### Step (6) : get genotypes for drought panel

2,0 for homs reference
1,1 for hets
0,2 for homs alternative

[vcf_to_read_counts.awk](https://github.com/rozenn-pineau/drought-project/blob/main/vcf_to_read_counts.awk)

### Step (7) : put the file together

Paste columns together to make the full file, following instructions on ancestry_hmm github page. 



### Processing results from ancestry_hmm output

Ancestry_hmm gives as an output one file per sample, that has the following header :
```
chrom	position	2,0,0	1,1,0	1,0,1	0,2,0	0,1,1	0,0,2
```
"2,0,0" column has the probability that the loci was homozygote for var. tuberculatus before hybridization, and stayed homozygote for var. tuberculatus at pulse 2 and 3 (0 chromosomes from var. rudis). GT = 0

"1,1,0" column has the probability that the loci was heterozygote for var. tuberculatus before hybridization, and became heterozygote for var. rudis at pulse 2 (1 chromosome from var. rudis). GT = 1

"1,0,1" column has the probability that the loci was heterozygote for var. tuberculatus before hybridization, and became heterozygote for var. rudis at pulse 3 (1 chromosome from var. rudis). GT = 1

"0,2,0" column has the probability that the loci was homozygote for var. rudis before hybridization, and stayed homozygote for var. rudis (0 chromosome from var. tuberculatus). GT = 2

"0,1,1" column has the probability that the loci was homozygote for var. tuberculatus before hybridization, and became gain one var. rudis loci at each pulse (1 chromosome from var. rudis). GT = 2

"0,0,2" column has the probability that the loci was homozygote for var. rudis before hybridization, and became homzygote for var. rudis at pulse 3 (1 chromosome from var. rudis). GT = 2


The probabilitiees in the columns add up to 1. 

We choose a threshold of 0.8 to assign ancestry at each site using the following script. 


### Step (1) : gather all samples in one file

```
#!/bin/bash

folder_name=$(basename "$PWD")
output_file="${folder_name}_values.txt"

#initialize dataset with chrom and pos
awk 'BEGIN { OFS="\t" } {print $1,$2}' P11_Ag_11_T.posterior > $output_file

#loop through each individual and paste information to previous version of the file
for file in *.posterior; do

    header_name=$(basename "$file" .posterior)
    
    awk -v header="$header_name" 'NR == 1 { print header; }
    
    NR > 1 {
    
    result = "NA";
    if ($3 > 0.8) result = "0";
    else if ($4 > 0.8) result = "1";
    else if ($5 > 0.8) result = "1";
    else if ($6 > 0.8) result = "2";
    else if ($7 > 0.8) result = "2";
    else if ($8 > 0.8) result = "2";
    print result;
    }' $file > tmp
    
    paste $output_file tmp > tmp2
    
    mv tmp2 $output_file

done

#check: file is 786262 lines and 284 columns (chromosome, position and 282 individuals).

```

### Visualize : plot ancestry by chromosome in R

We used the following script [plot_ancestry.R](https://github.com/rozenn-pineau/drought-project/blob/main/plot_ancestry.R) to visualize ancestry per chromosome for each smaple in the dataset, ordered by longitude. 

To make sure independent runs gave consistent ancestry calls, I compared the output from two runs that had different starting proportions (-.3 and -.37 versus  -.5 and -.17) : [compare_ancestry_hmm_runs.R](https://github.com/rozenn-pineau/Drought-paper/blob/main/compare_ancestry_hmm_runs.R). I tested 1000 random sites on two chromosomes and dit not identify differences in ancestry calls between models, suggesting the robustness of the results. 


