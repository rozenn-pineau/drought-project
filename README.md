# drought-project
Scripts and analyses for the drought paper.


# Outline
[Figure 1 analyses](#Figure-1-analyses)

[Trade-offs analyses](#Trade-offs-analyses)

  [Herbicide resistance estimations](###Herbicide-resistance-estimations)


[Figure 2 analyses](#Figure-2-analyses)

[Figure 3 analyses](#Figure-3-analyses)

[Figure 4 analyses](#Figure-4-analyses)

[Figure 5 analyses](#Figure-5-analyses)

[Figure 6 analyses](#Figure-6-analyses)


# Figure 1 analyses
Drought LD50 calculations, multiple linear regression models and Wilcoxon tests : [figure1_v2.R](https://github.com/rozenn-pineau/drought-project/figure1_v2.R) and [Figure1E.Rmd](https://github.com/rozenn-pineau/drought-project/Figure1E.Rmd).

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
Rscript to plot results from trade-offs analyses : [trade_offs_control.Rmd](https://github.com/rozenn-pineau/drought-project/trade_offs_control.Rmd)

# Figure 3 analyses
Drought LD50 calculations, multiple linear regression models and Wilcoxon tests : [figure1_v2.R](https://github.com/rozenn-pineau/drought-project/edit/main/figure1_v2.R) and [Figure1E.Rmd](https://github.com/rozenn-pineau/drought-project/edit/main/Figure1E.Rmd).
