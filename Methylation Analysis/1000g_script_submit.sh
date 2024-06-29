#!/bin/bash

#Script for running the bedtools interect and quantification script on each 1000g sample

script=1000g_modkit_intersect.sh
RScript=calcMeth_new.R
#sample_processed=$1
WORKING_DIR=$1

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
