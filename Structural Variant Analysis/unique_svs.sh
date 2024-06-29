#!/bin/bash

# .txt file list of file paths to all preprocessed, unzipped Napu output sniffles2 vcfs in reverse order of how you want them in the plot (ie. AFR samples at right end of the plot = AFR samples first in the list)
file_paths_txt="list.txt"

# Make iterative files (samples to merge adding one sample at a time)
# Iteration_001 will be the last sample in list.txt, iteration_002 will be the last two samples, etc. 

mapfile -t filenames < list.txt
total_files=${#filenames[@]}

for ((i = 1; i <= total_files; i++)); do
    # Pad the index with leading zeros
    index_padded=$(printf "%03d" $i)
    new_filename="iteration_$index_padded.txt"

    lines_to_read=$i

    tail -n "$lines_to_read" list.txt > "$new_filename"
done

#---------#

##Run jasmine on all file iterations

for file in iteration_*; do
        number=$(echo "$file" | sed 's/iteration_//')
        jasmine_output="merge_$number"
        jasmine file_list=$file out_file=$jasmine_output --allow_intrasample --output_genotypes --ignore_strand --dup_to_ins --centroid_merging genome_file=hg38.no_alt.fa
done

#----------#

##Count SVTYPE in each jasmine output

SVTYPE_count="SVTYPE_count.txt"

echo "SVTYPE Count" > "$SVTYPE_count"

for file in merge*; do
        bcftools query -f '%SVTYPE\n' "$file" | sort | uniq -c | awk -v filename="$file" '{ print $0, filename }'>> "$SVTYPE_count"
done

##fix header

sed -i '1s/^.*$/Count SVTYPE File/' SVTYPE_count.txt