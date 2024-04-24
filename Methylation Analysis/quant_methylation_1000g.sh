#!/bin/bash

#Processes 1000g outputs for quant_methylation. Will only run on card directories that have not been previously processed

INPUT_DIR=$1
TEXT_DIR=$2 #directory for list of processed samples
OUTPUT_DIR=$3
REFERENCE=$4

module load modkit/0.1.11
module load  samtools/1.17

cd $TEXT_DIR

rm samples_already_processed.txt # remove old list to regenerate new list for each new run

log_file="quant_methylation_1000g_log.txt"

exec >> "$log_file" 2>&1

# Output the current list of samples processed. This is used for the intersect script to only apply to newly processed directories
output_file="samples_already_processed.txt"

# Loop through subdirectories and write them to the output file
for dir in "$OUTPUT_DIR"/*/; do
    basename "${dir%/}" >> "$output_file"
done

readarray -t subdirs < "$output_file"


cd $OUTPUT_DIR

# Loop through each directory in the input directory
for dir in $INPUT_DIR/*; do
    # Get the basename of the directory
    dirname=$(basename "${dir%/}")

    # Check if the basename is NOT in the list of output directory names
    if [[ ! " ${subdirs[@]} " =~ " $dirname " ]]; then
        # Run modkit to generate bedMethyl pileup files. Currently only analyzing 5mC
        file=$(find $INPUT_DIR/$dirname -type f -name "*.haplotagged.bam")
        modkit pileup "$file" --cpg --ref $REFERENCE -t 40 --ignore h --combine-strands --partition-tag HP "$OUTPUT_DIR/${dirname}" --log-filepath "$OUTPUT_DIR/${dirname}/pileup.log" --prefix "${dirname}.cpg"
        echo "sample complete"
    fi
done


echo "modkit complete"