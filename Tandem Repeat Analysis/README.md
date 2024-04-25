These are scripts for performing the tandem repeat analysis

**Part One Assembling the bed file for geneotyping**

See the shell script "" for how the loci genotyped were selected


**Part Two: Run Vamos**

Make sure vamos is installed:

The conditions vamos was run is in the shell script "". A for-loop iterates through the hapdiff assembly files generated from the Napu pipeline (). These files are on AWS ().


**Part Three: Extract genotype information**

The python script "" counts the motif number at each locus and combined the information into a single tab-separated value. Sample vcfs can be combined into a single vcf using the combine_vcf.py script which was written by Ren et al., for vamos. This format is more difficult to parse.

**Part Four: Analysis and plots**

The tsv can be analysed in R or Python.
