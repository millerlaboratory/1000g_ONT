#!/bin/bash
#Use bedtools to filter for CpGs found in CpG islands, add sample and haplotype labels to each file and calculate the mean fraction of reads at each CpG island on the X chromosome

WORKING_DIR=$1
OUTPUT_DIR=$2
Reference_bed=/n/users/sgibson/1000g_preprint/methylation/imprinting_loci_for_intersect.bed

cd $WORKING_DIR

module load bedtools/2.30.0

sample=$(basename $WORKING_DIR | cut -d'.' -f1)

# Loop through each BED file in the directory for zipped files
for bed_file in $WORKING_DIR/*cpg_[1-2].bed.gz; do
    # Perform bedtools intersect with the current BED file
    sample=$(echo "${bed_file##*/}" | cut -d'.' -f1) #Sample labeling compliant with 1000g samples
    hap=$(echo "${bed_file##*/}" | cut -d'_' -f2 | cut -d '.' -f1) #Haplotype scheme compliant with 1000g naming structure
    bedtools intersect -wa -wb -a "$bed_file" -b $Reference_bed \
    | awk -v OFS="\t" -v sample="$sample" -v hap="$hap" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\t" $20 "\t" $21 "\t" $22 "\t" $23 "\t" sample "\t" hap}' > "$OUTPUT_DIR/${sample}_${hap}_intersected.bed"
done

cd $OUTPUT_DIR

#Merge all of the samples into one bed file
cat *_intersected.bed > "${sample}_dmr.bed"

rm *_intersected.bed


echo "$WORKING_DIR processed"