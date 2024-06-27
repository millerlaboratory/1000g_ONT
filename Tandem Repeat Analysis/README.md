These are scripts for performing the tandem repeat analysis

**Part One Assembling the bed file for geneotyping**

See the shell script  for how the loci genotyped were selected, the original motif set can be found here: https://zenodo.org/records/835736


**Part Two: Run Vamos**

Make sure vamos is installed:

The conditions vamos was run is in the shell script **vamos_original.sh**. A for-loop iterates through the hapdiff assembly files generated from the Napu pipeline (https://github.com/nanoporegenomics/napu_wf).


**Part Three: Extract genotype information**

The python script **count.py** counts the motif number at each locus and combined the information into a single tab-separated value. Sample vcfs can be combined into a single vcf using the **combine_vcf.py** script which was written by Ren et al., for vamos.

An example bash wrapper can be found too (**count_wrapper.sh**)

Disease-associated loci can be filtered like the method in **filter_STR.sh** This bed file is in the supplemental tables.

**Part Four: Analysis and plots**

The tsv can be analysed in R or Python. Figures from the paper are in Fig5.R
