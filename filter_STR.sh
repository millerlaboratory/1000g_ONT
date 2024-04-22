#!/bin/bash

#For manuscript analysis, filter sample_allele vcfs for disease loci

VCF_DIR=/n/users/sgibson/Projects/1KGP-ONT/1000g_preprint/all_simple_repeats/vamos_sr_vcf_original
BED=/n/users/sgibson/Projects/1KGP-ONT/1000g_preprint/all_simple_repeats/STR-disease-loci.processed_revised.bed
OUT_DIR=/n/users/sgibson/Projects/1KGP-ONT/20240422_Revision_Scripts/filtered_vcf/

module load bedtools/2.30.0

for file in $VCF_DIR/*.vcf; do
    sample=$(basename $file | cut -d "_" -f 1)
    hp=$(basename $file | cut -d "_" -f 5 | cut -d "." -f 1)
    bedtools intersect -wa -a $file -b $BED > $OUT_DIR/${sample}_${hp}_STR-disease-loci.vcf
done