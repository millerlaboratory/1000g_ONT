#!/bin/bash

module load bcftools
module load bedtools

# Starting with the inter-sample jasmine merged intra-sample confident call sets per individual (20240423_jasmine_intrasample_noBND_custom_suppvec_alphanumeric_header_JASMINE.vcf.gz)
vcf="20240423_jasmine_intrasample_noBND_custom_suppvec_alphanumeric_header_JASMINE.vcf.gz"


# Parse jasmine_fixed_lengths.bed for SVs that intersect OMIM exons
OMIM_EXONS="OMIM_exons_hg38_sorted.bed"
bcftools view -R $OMIM_EXONS "$vcf" > jasmine_output_OMIM_exons.vcf

# Query the vcf for relevant fields 
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVLEN\t%INFO/SVTYPE\t%INFO/SUPP\n' jasmine_output_OMIM_exons.vcf > Figure_S12_data.txt