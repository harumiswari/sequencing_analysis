#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 17:19:54 2024

@author: gwisna
"""

import pandas as pd
import pysam

def categorize_all_reads(samfile, target_chrom_1, target_chrom_2):
    all_reads = []
    read_names_processed = set()

    def process_reads(chromosome):
        for read in samfile.fetch(chromosome):
            if read.query_name in read_names_processed or read.is_unmapped:
                continue
            read_names_processed.add(read.query_name)
            original_chrom = read.reference_name
            original_pos = read.reference_start
            cigar = read.cigarstring
            read_details = (read.query_name, f"{original_chrom}:{original_pos}", cigar)
            all_reads.append(read_details)

    process_reads(target_chrom_1)
    process_reads(target_chrom_2)

    return all_reads

def process_barcode(barcode):
    base_path = f'/Users/gwisna/Desktop/hF9mousecombined/BC{barcode}/'
    bam_path = f'{base_path}BC{barcode}_deduplicated.bam'
    samfile = pysam.AlignmentFile(bam_path, "rb")

    target_chrom_1 = 'CM001000.3'
    target_chrom_2 = 'ITRtoITRinsert'

    all_reads = categorize_all_reads(samfile, target_chrom_1, target_chrom_2)

    # Convert to DataFrame and deduplicate
    all_reads_df = pd.DataFrame(all_reads, columns=['Read Name', 'Original Align', 'Cigar'])
    all_reads_df.drop_duplicates(subset='Read Name', keep='first', inplace=True)

    # Extract chromosome and position
    all_reads_df['Chrom'] = all_reads_df['Original Align'].apply(lambda x: x.split(':')[0])
    all_reads_df['Position'] = all_reads_df['Original Align'].apply(lambda x: int(x.split(':')[1].split(' ')[0]))

    # Process primary data
    primary_data = all_reads_df[all_reads_df['Chrom'] == target_chrom_1]
    primary_location_counts = primary_data.groupby(['Chrom', 'Position']).size().reset_index(name='Counts')
    primary_output_path_counts = f'{base_path}BC{barcode}_{target_chrom_1}_reads_counts.csv'
    primary_location_counts.to_csv(primary_output_path_counts, index=False)
    primary_data[['Read Name']].to_csv(f'{base_path}BC{barcode}_{target_chrom_1}_read_names.txt', index=False, header=False)

    # Process secondary data
    secondary_data = all_reads_df[all_reads_df['Chrom'] == target_chrom_2]
    secondary_location_counts = secondary_data.groupby(['Chrom', 'Position']).size().reset_index(name='Counts')
    secondary_output_path_counts = f'{base_path}BC{barcode}_{target_chrom_2}_reads_counts.csv'
    secondary_location_counts.to_csv(secondary_output_path_counts, index=False)
    secondary_data[['Read Name']].to_csv(f'{base_path}BC{barcode}_{target_chrom_2}_read_names.txt', index=False, header=False)

    samfile.close()

    print(f"Read name lists and counts from {target_chrom_1} saved to {primary_output_path_counts}")
    print(f"Read name lists and counts from {target_chrom_2} saved to {secondary_output_path_counts}")

# Loop over barcodes
for barcode in range(11, 21):
    process_barcode(barcode)


######To get fastq file from filtered bam#######

def filter_fastq(fastq_file, read_names, output_file):
    """Filter a FASTQ file based on a set of read names and write the matching records to a new file."""
    with open(fastq_file, 'r') as fq, open(output_file, 'w') as out_fq:
        while True:
            header = fq.readline() 
            if not header:
                break 
            sequence = fq.readline() 
            plus = fq.readline()  
            quality = fq.readline() 

            read_name_in_header = header.split()[0][1:]  
            if read_name_in_header in read_names:
                out_fq.write(header + sequence + plus + quality)

# Loop over barcodes
for barcode in range(11, 21):
    data_path = f'/Users/gwisna/Desktop/hF9mousecombined/BC{barcode}/BC{barcode}_ITRtoITRinsert_read_names.txt'
    data = pd.read_csv(data_path, delimiter='\t', header=None, names=['Read Name'])
    read_names = set(data['Read Name'])

    fastq_file = f'/Users/gwisna/Desktop/with__pear/BC{barcode}/BC{barcode}_umi.fastq'
    output_file = f'/Users/gwisna/Desktop/hF9mousecombined/BC{barcode}/BC{barcode}_ITRtoITRinsert.fastq'

    filter_fastq(fastq_file, read_names, output_file)
