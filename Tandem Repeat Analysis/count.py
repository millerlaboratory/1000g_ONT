import os
import re
import pandas as pd
from collections import defaultdict 
import argparse

#Count the number of motifs in the sequence, and calculate the total length of the sequence and filter out second alleles for 46, XY indviduals on the X chromosome. The sequence annotation is included for making waterfall plots.

def calculate_total_length(sequence, vntr_to_motif_length):
    total_length = 0
    for motif_index in sequence:
        total_length += vntr_to_motif_length[motif_index]
    return total_length

def process_vcf_files(directory_path, bed_file, sample_metadata, output_path):
    vcf_files=[file for file in os.listdir(f"{directory_path}") if file.endswith(".vcf")]

    # Parse the motifs in the bed file and calculate their sizes
    vntr_to_motif_length = {}
    motif_number = 0

    bed=f"{bed_file}"
    with open(bed, 'r') as file:
        for line in file:
            data = line.strip().split('\t')
            chr = data[0]
            start = data[1]
            end = data[2]
            vntr_bed = chr + ":" + start + "-" + end
            motifs = data[3].split(',')
            vntr_to_motif_length[vntr_bed] = [len(motif) for motif in motifs]

    results = []
    for vcf in vcf_files:
        file = f"{directory_path}/{vcf}"
        with open(file, 'r') as f:
            lines = f.readlines()
            for idx, line in enumerate(lines):
                if line.startswith("#"):
                    continue
                fields = re.split('\t', line.strip('\n'))
                sample_name = vcf.split("_")[0]
                allele=vcf.split("_")[1]
                chrm = fields[0]
                start = fields[1]
                infos = re.split(";", fields[7])
                end = re.split("=", infos[0])[1]
                vntr = chrm + ":" + start + "-" + end
                anno = re.split("=", infos[3])[1]
                sequence = list(map(int, anno.split(',')))
                count=len(sequence)

                # Calculate the length for the current line
                vntr_motif_length = vntr_to_motif_length.get(vntr, [0 for _ in sequence])
                total_length = calculate_total_length(sequence, vntr_motif_length)
                results.append({'chr' : chrm, 'start' : start, 'end' : end, 'sample' : sample_name, 'vntr': vntr, 'count': count, 'length': total_length, 'anno': anno, 'allele': allele})

    # Create a dataframe and write the results
    df = pd.DataFrame(results)

    sample_metadata = pd.read_csv(sample_metadata, sep='\t', header=None, names=['sample', 'sex', 'population'])

    df = df.merge(sample_metadata, on='sample')

    df = df[~((df['sex'] == 'male') & (df['chr'] == 'chrX') & (df['allele'] == 'hp2'))]

    df.to_csv(f'{output_path}/first100_strchive_count.tsv', sep='\t', index=False)

    return df





    

def main():
    parser = argparse.ArgumentParser(description='Process allele data.')
    parser.add_argument('input_path', help='Path to the input data')
    parser.add_argument('output_path', help='Path to the output CSV file')
    parser.add_argument('bed_file', help='bed file with motifs')
    parser.add_argument('sample_metadata', help='Metadata on the samples')
    args = parser.parse_args()

    df = process_vcf_files(args.input_path, args.bed_file, args.sample_metadata, args.output_path)


if __name__ == "__main__":
    main()

