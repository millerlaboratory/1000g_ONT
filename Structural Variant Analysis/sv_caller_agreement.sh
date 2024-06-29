#!/bin/bash

# Directory containing all SV vcfs for all individuals:
main_directory="/all/vcfs"

# Move files into subdirectories based on their Sample ID
for file in "$main_directory"/*vcf*; do
    if [ -f "$file" ]; then
        # Get the first 7 characters of the file name
        prefix=$(basename "$file" | cut -c 1-7)

        subdirectory="$main_directory/$prefix"

        # Check if the subdirectory exists, if not, create it
        if [ ! -d "$subdirectory" ]; then
            mkdir -p "$subdirectory"
        fi

        # Move the file into the corresponding subdirectory
        mv "$file" "$subdirectory/"
    fi
done

# Create the file path lists for all vcfs per individual 

for subdir in "$main_directory"/*/; do
        find "$subdir" -type f -name "*.vcf" > "$subdir"/list.txt
done

# Run Jasmine merging for all vcfs per individual

for subdir in "$main_directory"/*/; do
        subdirname=$(basename "$subdir")
        jasmine file_list="$subdir"/list.txt out_file="$subdir"/"$subdirname"_jasmine_intrasample.vcf threads=50 --allow_intrasample --output_genotypes --ignore_strand --dup_to_ins --centroid_merging genome_file=hg38.no_alt.fa
done

# Move all jasmine output files to a new directory 

find $main_directory -type f -name "*jasmine_intrasample.vcf" -exec scp {} /jasmine/output/ \;

# query the jasmine output vcfs for relevant info --> .txt file and concatenate 

for vcf_file in "$vcf_directory"/*.vcf; do
    bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/END\t%INFO/SVLEN\t%INFO/SVTYPE\t%INFO/SUPP_VEC\n' "$vcf_file" > "${vcf_file%.vcf}.txt"
done

cat /jasmine/output/*.txt > Figure_3E_S9_data.txt

# Remove BNDs from jasmine merge outputs

for vcf_file in /jasmine/output/*vcf; do
    awk '!/BND/' "$vcf_file" > "${vcf_file%.vcf}_noBND.vcf"
done


# Parse jasmine output vcfs (without BNDs) for SUPP_VECs that represent "confident" combinations of callers (for the maniuscript, we used any SV that was called by hapdiff and at least 2 unique alignment-based callers)

for file in /jasmine/output/*noBND.vcf; do
        output_file="${file%.vcf}_custom_suppvec"
        bcftools view -i 'INFO/SUPP_VEC="11111" || INFO/SUPP_VEC="11110" || INFO/SUPP_VEC="11101" || INFO/SUPP_VEC="11011" || INFO/SUPP_VEC="11100" || INFO/SUPP_VEC="11001" || INFO/SUPP_VEC="10111"|| INFO/SUPP_VEC="10110" || INFO/SUPP_VEC="10101" | INFO/SUPP_VEC="10011"' -o "$output_file" "$file"
done

## Run intersample merge on intrasample confident callsets 
# Make a .txt list of all of the file paths to the confident vcfs
find /jasmine/output/*_custom_suppvec.vcf > custom_suppvec_list.txt

# Run Jasmine
jasmine file_list=custom_suppvec_list.txt out_file=custom_suppvec_merged.vcf threads=50 --allow_intrasample --output_genotypes --ignore_strand --dup_to_ins --centroid_merging genome_file=hg38.no_alt.fa

# Query the vcf
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVLEN\t%INFO/SVTYPE\t%INFO/SUPP\t%INFO/IDLIST\n' custom_suppvec_merged.vcf > custom_suppvec_merged.txt 

# Parse the output for variants seen in exactly 2 samples
awk '$6 == 2' custom_suppvec_merged.txt > custom_suppvec_merged_exactly2samples.txt

# Split the IDLIST column by comma
awk -F'\t' '{
  split($7, a, ",");
  print $1, $2, $3, $4, $5, $6, a[1], a[2], $8, $9, $10
}' OFS='\t' custom_suppvec_merged_exactly2samples.txt > custom_suppvec_merged_exactly2samples_2.txt

# Parse the output for the ancestry of the 2 samples (characters 12-14 of the sample ID)
awk -F'\t' '{$7 = substr($7, 12, 3); $8 = substr($8, 12, 3); print}' OFS='\t' custom_suppvec_merged_exactly2samples_2.txt > custom_suppvec_merged_exactly2samples_ANC.txt

# Combine ancestry info in alphabetical order 
awk -F'\t' '{if ($7 < $8) {$7 = $7 "_" $8;} else {$7 = $8 "_" $7;} $8 = ""; print}' OFS='\t' custom_suppvec_merged_exactly2samples_ANC.txt > Figure_3F_data.txt
