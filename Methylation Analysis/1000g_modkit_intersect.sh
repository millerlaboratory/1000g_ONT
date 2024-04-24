#!/bin/bash
#Use bedtools to filter for CpGs found in CpG islands, add sample and haplotype labels to each file and calculate the mean fraction of reads at each CpG island on the X chromosome

WORKING_DIR=$1
#RScript=$2
Reference_bed=/n/users/sgibson/reference/chrX_cpgIslands.hg38_NUM_FOR_INTERSECT.bed #Some error submitting from the script submit

cd $WORKING_DIR

module load bedtools/2.30.0
module load R/4.2.3

# Loop through each BED file in the directory for zipped files
for bed_file in $WORKING_DIR/*cpg_[1-2].bed.gz; do
    # Perform bedtools intersect with the current BED file
    bedtools intersect -wa -wb -a "$bed_file" -b $Reference_bed \
    | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $22}' > "${bed_file%.bed.gz}_Islands.bed"
done

#Add the sample name and haplotype as columns in each file
for bed_file in $WORKING_DIR/*_Islands.bed; do
    sample=$(echo "${bed_file##*/}" | cut -d'.' -f1) #Sample labeling compliant with 1000g samples
    hap=$(echo "${bed_file##*/}" | cut -d'_' -f2) #Haplotype scheme compliant with 1000g naming structure
    awk -v OFS="\t" -v sample="$sample" -v hap="$hap" '{print $0 "\t" sample "\t" hap}' "$bed_file" > "${bed_file%.bed}_labeled.bed"
done

#Merge all of the samples into one bed file
cat *_labeled.bed > chrX.bed

rm *_labeled.bed
rm *_Islands.bed
rm merged_haplotypes.bed

#Split the combined, processed file by chromosome
#awk '{print > $1".bed"}' merged_haplotypes.bed

#Runs the R script to calculate the mean fraction of methylated reads for each CpG island on the X chromosome
#Rscript $RScript $WORKING_DIR chrX 

#Compress whole-genome bed files
#bed_files=$(find $INPUT_DIR -type f -name "*cpg_[1-2].bed")

#for file in $bed_files; do
  #bgzip "$file"
#done

#bgzip *_ungrouped.bed

echo "$WORKING_DIR processed"