#!/bin/bash

BAM1=HG02389.LSK110.R9_PMDV_FINAL.haplotagged.bam
BAM2=HG03022.LSK110.R9_PMDV_FINAL.haplotagged.bam
GTF=gencode.v38.annotation.sorted.gtf.gz
OUT=$1

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
