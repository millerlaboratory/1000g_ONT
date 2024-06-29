#!/bin/bash

module load bcftools/1.19 #https://samtools.github.io/bcftools/howtos/index.html
module load hap.py/0.3.14 #https://github.com/Illumina/hap.py

# download merged llumina vcf for chr21 (https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz)

#List all of the samples in the illumina vcf
bcftools query -l illumina_chr21.vcf.gz > illumina_chr21_sample_list.txt

# -c only outputs variants seen in the specific sample (otherwise you get outputs for all variants seen in all samples in the jopint call, with 0/0 genotypes for variants not seen in specific samples - makes vcfs way larger)
# All samples from the 1000G Illumina dataset are already filtered for FILTER = PASS, so there is no "PASS" flag to match the filtering for the query files
while IFS= read -r sample; do
        output_file="${sample}.illumina_chr21_c1_PASS.vcf.gz"
        bcftools view -s "${sample}" illumina_chr21.vcf.gz -c 1 -Oz -o "${output_file}"
        tabix -p vcf "${output_file}"
done < illumina_chr21_sample_list.txt


# query directory is the directory containing the ONT clair3 and PMDV vcfs
query_directory=""

# truth directory is the dirctory containing the sample-parsed illumina vcfs 
truth_directory=""

# output directory is where the hap.py output goes
output_directory=""

# GIAB-defined high confidence regions bed file (https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed)
high_conf="HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

# hg38 FASTA file as suggested by https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use
hg38_fasta="hg38.no_alt.fa"

for subdir in "$query_directory"/; do
# Loop through the files in the query directory
        for query_file in "$subdir"/*.vcf.gz; do
        # Extract the prefix of the query file
        query_prefix=$(basename "$query_file" | cut -c 1-7)
        output_file_long=$(basename "$subdir")

        output_file="$output_directory/${output_file_long}"        
        truth_file=$(find "$truth_directory" -name "$query_prefix*gz" -type f)
              
        # Run hap.py command
        hap.py "$truth_file" "$query_file" -T "$high_conf" -o "$output_file"  -r "$hg38_fasta"
done