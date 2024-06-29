#!/bin/bash

#Script with bedtools commands to create the motif reference bedfile used for STR/VNTR genotyping

module load bedtools/2.30.0

#Filter the original motif bed from Ren et al., 2023 for the high confience regions
bedtools intersect -wa -a q-0.1_motifs.set148.bed -b FROM_MISHA_GIAB_SV_Tier1_grch38.fixed.bed > original_motifs.set148_GIAB_SV_Tier1_liftover_grch38.bed

#Grab the simple repeats from the repeat track from the UCSC genome browser
awk '$4 ~ /Simple_repeat/' repeats.hg38.sorted.bed > simple_repeats.hg38.sorted.bed

#Filter the motif bed for simple repeats
bedtools intersect -wa -a original_motifs.set148_GIAB_SV_Tier1_liftover_grch38.bed -b simple_repeats.hg38.sorted.bed > original_motifs.set148_GIAB_SV_Tier1_liftover_grch38_simple.bed

#Repace the STRchive loci with the specified motif information
bedtools intersect -v -a original_motifs.set148_GIAB_SV_Tier1_liftover_grch38_simple.bed -b STR_disease_loci_for_vamos.bed > original_motifs.set148_GIAB_SV_Tier1_liftover_grch38_simple_noSTRchive.bed


cat original_motifs.set148_GIAB_SV_Tier1_liftover_grch38_simple_noSTRchive.bed STR-disease-loci.processed_revised.bed > original_motifs.set148_GIAB_SV_Tier1_liftover_grch38_simple_STRchive.bed


bedtools sort -i original_motifs.set148_GIAB_SV_Tier1_liftover_grch38_simple_STRchive.bed > original_motifs.set148_GIAB_SV_Tier1_liftover_grch38_simple_STRchive_sorted_ARfix.bed
