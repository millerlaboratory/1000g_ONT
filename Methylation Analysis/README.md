Scripts and supplemental bed files for performing analysis

**Part One: Haplotype-Resolved Methylation Pileups**
Run modkit, specific parameters in **quant_methylation_1000g.sh**


**Part Two: X Chromosome Analysis**
Intersect for CpGs on the X-chromosome and calculate the average at each CpG island with **1000g_modkit_intersecy.sh** and **calcMeth_new.R**. These are run per-sample with the **1000g_script_submit.sh** wrapper

Further X-chromsome analysis is in **Fig6.R**

**Part Three: Phased Autosomal DMR analysis**
Filter for the autosome CpGs with the **get_pileup.sh** script and calculate averages with **mean_per_dmr.R**. This process is nearly identical to the chrX analysis, only different in when it combines the samples into one file.

The analysis is in **Fig6.R**


**Part Four: MeOW and Downstream Analysis**

MeOW: https://github.com/mgaleyuw/MeOW

PCA using the pileup was performed in the **PCA.R** script
