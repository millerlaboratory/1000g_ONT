#!/bin/bash

#Filtered SV vcfs for variants that are >=50bp, PASS quality filters, and are on full chromosomes 1-22,X,Y,M
  
for vcf in *sniffles2*vcf.gz; do
bcftools view -i '(INFO/SVTYPE="BND") || (INFO/SVTYPE="INS" || INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP" || INFO/SVTYPE="INV") && (INFO/SVLEN > 49 || INFO/SVLEN < -49)' -f PASS -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM -o PASS_chr_len_"$vcf" "$vcf"
done

for vcf in *cuteSV*vcf.gz; do
bcftools view -i '(INFO/SVTYPE="BND") || (INFO/SVTYPE="INS" || INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP" || INFO/SVTYPE="INV") && (INFO/SVLEN > 49 || INFO/SVLEN < -49)' -f PASS -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM -o PASS_chr_len_"$vcf" "$vcf"
done

for vcf in *hapdiff*vcf.gz; do
bcftools view -i '(INFO/SVTYPE="BND" || INFO/SVTYPE="INV") || (INFO/SVTYPE="INS" || INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP") && (INFO/SVLEN > 49 || INFO/SVLEN < -49)' -f PASS -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM -o PASS_chr_len_"$vcf" "$vcf"
done

for vcf in *svim*vcf.gz; do
sed -i 's/DUP:INT/DUP/g' \ sed -i 's/DUP:TANDEM/DUP/g' \ bcftools view -i 'QUAL > 10' \ bcftools view -i '(INFO/SVTYPE="BND" || INFO/SVTYPE="INV") || (INFO/SVTYPE="INS" || INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP") && (INFO/SVLEN > 49 || INFO/SVLEN < -49)' -f PASS -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM -o PASS_chr_len_"$vcf" "$vcf"
done