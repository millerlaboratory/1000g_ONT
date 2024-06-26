#loading flye 
module load flye/2.9.2
module load minimap2/2.26

# Loop through all FASTA files in the current directory
for fasta_file in *.fasta; do
    # Check if the file exists and is not a directory
    if [ -f "$fasta_file" ]; then
        # Extract the filename without extension
        file_name_without_ext="${fasta_file%.fasta}"
        
        # Run flye with the specified options
        flye --threads 50 --nano-hq "$fasta_file" --out-dir "$file_name_without_ext"
        
        # Optional: You can add further processing or actions here if needed
    fi
done