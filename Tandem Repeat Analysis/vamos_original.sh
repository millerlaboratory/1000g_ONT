#!/bin/bash

INPUT_DIR=/n/1000g/align-card-2.24-hg38/FIRST_100/
OUTPUT=/n/users/sgibson/Projects/1KGP-ONT/20240422_Revision/REVISED_S2/vamos_original_vcf
MOTIFS=/n/users/sgibson/Projects/1KGP-ONT/20240422_Revision/REVISED_S2/original_motifs.set148_GIAB_SV_Tier1_liftover_grch38_simple_STRchive_sorted_ARfix.bed

cd $OUTPUT


#MODULES
module load vamos/1.3.6
module load python/3.9

#VAMOS
for dir in $INPUT_DIR/*; do
    # Get the basename of the directory
    dirname=$(basename "${dir%/}")
    filename=$(echo $dirname | cut -d '-' -f1)
    hp1=$(find $INPUT_DIR/$dirname -type f -name "*.hapdiff_mat.bam")
    hp2=$(find $INPUT_DIR/$dirname -type f -name "*.hapdiff_pat.bam")
    vamos --contig -b $hp1 -r $MOTIFS -s $filename -o "$OUTPUT/${filename}_hp1.vcf" -t 8
    vamos --contig -b $hp2 -r $MOTIFS -s $filename -o "$OUTPUT/${filename}_hp2.vcf" -t 8
    echo "sample complete"
done
echo "Vamos complete"
