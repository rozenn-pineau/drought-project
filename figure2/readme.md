# Figure 2 analyses

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

Rscript to plot results admixture : [plot_PCA_res.R](https://github.com/rozenn-pineau/drought-project/blob/main/figure2/plot_PCA_res.R)
