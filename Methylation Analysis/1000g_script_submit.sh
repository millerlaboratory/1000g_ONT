#!/bin/bash

#Script for running the bedtools interect and quantification script on each 1000g sample

script=/n/users/sgibson/github_repo/quant_methylation/1000G_ONT_WORKFLOW/1000g_modkit_intersect.sh
RScript=/n/users/sgibson/github_repo/quant_methylation/1000G_ONT_WORKFLOW/calcMeth_new.R
#sample_processed=$1
WORKING_DIR=/n/users/sgibson/1000g_methylation/modkit_v0.1.11/FIRST_100_MODKIT

module load R/4.2.3

cd $WORKING_DIR

#readarray -t sample_done < "$sample_processed"

# Match the basename for each modkit processed directory with the list of ones that had already been processed
for subdir in "$WORKING_DIR"/*; do
    if [ -d "$subdir" ]; then
        dirname=$(basename "$subdir")
            cd "$subdir"
            Rscript $RScript $subdir chrX
            echo "Script executed for directory: $subdir"
    fi
done
