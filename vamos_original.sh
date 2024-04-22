#!/bin/bash
INPUT_DIR=/n/1000g/align-card-2.24-hg38/FIRST_100/
OUTPUT=/n/users/sgibson/1000g_preprint/all_simple_repeats/vamos_sr_vcf_original
MOTIFS=/n/users/sgibson/1000g_preprint/all_simple_repeats/original/original_motifs.set148_GIAB_SV_Tier1_liftover_grch38_simple_STRchive_sorted.bed

module load vamos/1.3.6
module load python/3.9

cd $OUTPUT

# Loop through each directory in the input directory
for dir in $INPUT_DIR/*; do
    # Get the basename of the directory
    dirname=$(basename "${dir%/}")
    filename=$(echo $dirname | cut -d '.' -f1)
    hp1=$(find $INPUT_DIR/$dirname -type f -name "*_hapdiff_mat.bam")
    hp2=$(find $INPUT_DIR/$dirname -type f -name "*_hapdiff_pat.bam")
    vamos --contig -b $hp1 -r $MOTIFS -s $filename -o "$OUTPUT/${filename}_GIAB_simple_original_hp1.vcf" -t 8
    vamos --contig -b $hp2 -r $MOTIFS -s $filename -o "$OUTPUT/${filename}_GIAB_simple_original_hp2.vcf" -t 8
    echo "sample complete"
done