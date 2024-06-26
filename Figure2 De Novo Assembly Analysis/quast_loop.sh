#!/bin/bash

#INPUT_DIR=
INPUT_DIR=
WORKING_DIR=
REF=GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta


module load quast

cd $WORKING_DIR

# Find all fasta files within subdirectories
fasta_files=$(find "$INPUT_DIR" -type f -name "assembly.fasta")

# Iterate over each fasta file
for file in $fasta_files; do
	file_name=$(basename"$file")
	dir_name=$(dirname "$file")             # Get the directory name
    folder_name=$(basename "$dir_name")      # Get the base name of the directory
    output_folder="$WORKING_DIR/$folder_name"
    quast.py "$file" \
	-r $REF \
	-t 10 \
	--large \
	--x-for-Nx 75 \
	--silent \
	--fragmented \
	-o "$output_folder"
done

