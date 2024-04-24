#!/bin/bash

BAM1=/n/1000g/align-card-2.24-hg38/FIRST_100/HG02389.LSK110.R9/HG02389.LSK110.R9_PMDV_FINAL.haplotagged.bam
BAM2=/n/1000g/align-card-2.24-hg38/FIRST_100/HG03022.LSK110.R9/HG03022.LSK110.R9_PMDV_FINAL.haplotagged.bam
GTF=/n/users/sgibson/1000g_preprint/methylation/figs/gencode.v38.annotation.sorted.gtf.gz
OUT=/n/users/sgibson/1000g_preprint/methylation/figs

cd $OUT

source activate test

modbamtools plot -r chr10:71372360-71382450 \
    --gtf $GTF \
    --track-titles Genes \
    --out $OUT \
	--samples HG02389,HG03022 \
    --hap \
    --prefix SLC29A3_plot_revised \
  $BAM1 $BAM2
