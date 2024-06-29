#!/bin/bash

VCF_DIR=$1
BED=STR-disease-loci.processed_revised.bed
OUT_DIR=$3
SAMPLE_METADATA=first100_1KGP_metadata.tsv #This is supplemental table 1 in the manuscript

cd $VCF_DIR

source activate jupyter

python3.12 count.py $VCF_DIR $OUT_DIR $BED $SAMPLE_METADATA
