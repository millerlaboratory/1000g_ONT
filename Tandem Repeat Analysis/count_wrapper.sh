#!/bin/bash

VCF_DIR=/n/users/sgibson/Projects/1KGP-ONT/20240422_Revision_Scripts/filtered_vcf
BED=/n/users/sgibson/Projects/1KGP-ONT/1000g_preprint/all_simple_repeats/STR-disease-loci.processed_revised.bed
OUT_DIR=/n/users/sgibson/Projects/1KGP-ONT/20240422_Revision_Scripts/count
SAMPLE_METADATA=/n/users/sgibson/Projects/1KGP-ONT/POST_100/repeat_outlier_tool/data/first100_1KGP_metadata.tsv #This is supplemental table 1 in the manuscript

cd $VCF_DIR

source activate jupyter

python3.12 /n/users/sgibson/Projects/1KGP-ONT/20240422_Revision_Scripts/count.py $VCF_DIR $OUT_DIR $BED $SAMPLE_METADATA
