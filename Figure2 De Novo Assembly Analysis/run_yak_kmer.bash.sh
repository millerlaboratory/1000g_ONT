#!/bin/bash
module load yak/0.1-r56


# Find all pairs of files matching the pattern
file_pairs=($(ls *_1.fastq.gz | sed 's/_1.fastq.gz//'))

# Loop through each pair of files
for name in "${file_pairs[@]}"; do
    file1="${name}_1.fastq.gz"
    file2="${name}_2.fastq.gz"
  
    
    # Run the yak command
    yak count -b37 -t50 -o "$name.yak" <(zcat < "$file1") <(zcat < "$file2")
done
