#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 18:14:31 2024

@author: harumi
"""

import pysam
import pandas as pd
from Bio import SeqIO

def categorize_chimeric_reads(samfile, chrom, valid_read_names):
    filtered_reads = []
    ITRtoITRinsert_reads = []
    read_names_processed = set()

    for read in samfile.fetch(chrom):
        if read.query_name in read_names_processed:
            continue
        if read.query_name not in valid_read_names:
            continue

        read_names_processed.add(read.query_name)
        original_chrom = read.reference_name
        original_pos = read.reference_start

        if read.has_tag('SA'):
            chimeric_info = read.get_tag('SA')
            for chimeric_part in chimeric_info.split(';'):
                if chimeric_part:
                    c_chrom, c_pos, c_strand, c_cigar, c_mapq, c_nm = chimeric_part.split(',')
                    if c_chrom != original_chrom:
                        filtered_reads.append((read.query_name, f"{original_chrom}:{original_pos}", f"{c_chrom}:{c_pos} {c_strand}", c_cigar))
                    if c_chrom == chrom:
                        ITRtoITRinsert_reads.append((read.query_name, f"{original_chrom}:{original_pos}", f"{c_chrom}:{c_pos} {c_strand}", c_cigar))
        else:
            if original_chrom == chrom:
                ITRtoITRinsert_reads.append((read.query_name, f"{original_chrom}:{original_pos}", f"{original_chrom}:{original_pos} {read.is_reverse}", read.cigarstring))
                
    return filtered_reads, ITRtoITRinsert_reads

def filter_fastq_and_extract_read_names(fastq_file, read_names):
    valid_read_names = set()
    target_sequence = "CAAGGGAACCTTGAGAGAG"
    excluded_sequence = "GGGGGGGGGG"

    with open(fastq_file, 'r') as fq:
        while True:
            header = fq.readline()
            if not header:
                break
            sequence = fq.readline().strip().upper()
            plus = fq.readline()
            quality = fq.readline()

            read_name_in_header = header.split()[0][1:]
            if read_name_in_header in read_names and excluded_sequence not in sequence and target_sequence in sequence:
                valid_read_names.add(read_name_in_header)
    return valid_read_names

def process_barcode(barcode):
    base_path = f'/Users/gwisna/Desktop/hF9mousecombined/BC{barcode}/'
    bam_path = f'{base_path}BC{barcode}_deduplicated.bam'
    samfile = pysam.AlignmentFile(bam_path, "rb")

    chromosome_name = 'ITRtoITRinsert'

    fastq_file = f'/Users/gwisna/Desktop/with__pear/BC{barcode}/BC{barcode}_umi.fastq'
    read_names_in_txt = set(pd.read_csv(f'{base_path}BC{barcode}_chimeric_reads_results_filtered.txt', delimiter='\t')['Read Name'])
    valid_read_names = filter_fastq_and_extract_read_names(fastq_file, read_names_in_txt)

    filtered_chimeric_reads, ITRtoITRinsert_reads = categorize_chimeric_reads(samfile, chromosome_name, valid_read_names)

    with open(f'{base_path}BC{barcode}_chimeric_reads_results_filtered.txt', 'w') as file:
        file.write("Read Name\tOriginal Align\tChimeric Align\tCigar\n")
        for read in filtered_chimeric_reads:
            file.write(f"{read[0]}\t{read[1]}\t{read[2]}\t{read[3]}\n")

    samfile.close()

    data_filtered = pd.read_csv(f'{base_path}BC{barcode}_chimeric_reads_results_filtered.txt', sep='\t')
    data_filtered.drop_duplicates(subset='Read Name', keep='first', inplace=True)
    data_filtered['Chrom'] = data_filtered['Chimeric Align'].apply(lambda x: x.split(':')[0])
    data_filtered['Position'] = data_filtered['Chimeric Align'].apply(lambda x: int(x.split(':')[1].split(' ')[0]))
    location_counts_filtered = data_filtered.groupby(['Chrom', 'Position']).size().reset_index(name='Counts')
    location_counts_filtered.to_csv(f'{base_path}BC{barcode}_integration_site_distribution_filtered.csv', index=False)

    output_fastq = f'{base_path}BC{barcode}_chimeric_reads_results_filtered.fastq'
    filter_fastq(fastq_file, valid_read_names, output_fastq)

def filter_fastq(fastq_file, valid_read_names, output_file):
    with open(fastq_file, 'r') as fq, open(output_file, 'w') as out_fq:
        while True:
            header = fq.readline()
            if not header:
                break
            sequence = fq.readline().strip().upper()
            plus = fq.readline()
            quality = fq.readline()

            read_name_in_header = header.split()[0][1:]
            if read_name_in_header in valid_read_names:
                out_fq.write(header + sequence + plus + quality)

for barcode in range(31, 41):
    process_barcode(barcode)

print("Processing complete for all files.")


######To extract the ITRtoITRinsert##########
import pysam
import pandas as pd
from Bio import SeqIO

def categorize_chimeric_reads(samfile, chrom, start, end, valid_read_names):
    filtered_reads = []
    ITRtoITRinsert_reads = []
    read_names_processed = set() 

    for read in samfile.fetch(chrom, start, end):
        if read.query_name in read_names_processed:
            continue 
        if read.query_name not in valid_read_names:
            continue 

        if read.has_tag('SA'):
            read_names_processed.add(read.query_name)
            original_chrom = read.reference_name
            original_pos = read.reference_start
            chimeric_info = read.get_tag('SA')
            for chimeric_part in chimeric_info.split(';'):
                if chimeric_part:
                    c_chrom, c_pos, c_strand, c_cigar, c_mapq, c_nm = chimeric_part.split(',')
                    if c_chrom != original_chrom:
                        filtered_reads.append((read.query_name, f"{original_chrom}:{original_pos}", f"{c_chrom}:{c_pos} {c_strand}", c_cigar))
                    if c_chrom == chrom:
                        ITRtoITRinsert_reads.append((read.query_name, f"{original_chrom}:{original_pos}", f"{c_chrom}:{c_pos} {c_strand}", c_cigar))
    return filtered_reads, ITRtoITRinsert_reads

def filter_fastq_and_extract_read_names(fastq_file, target_sequence, excluded_sequence):
    """Filter a FASTQ file based on the presence of target_sequence and absence of excluded_sequence. Return the valid read names."""
    valid_read_names = set()

    with open(fastq_file, 'r') as fq:
        while True:
            header = fq.readline()  
            if not header:
                break  
            sequence = fq.readline().strip().upper()  
            plus = fq.readline()  
            quality = fq.readline() 

            # Extract the read name from the header
            read_name_in_header = header.split()[0][1:]  # read name is the first part of the header, removes '@'
            if excluded_sequence not in sequence and target_sequence in sequence:
                valid_read_names.add(read_name_in_header)
                print(f"Valid read name: {read_name_in_header}")  
    return valid_read_names

def process_barcode(barcode):
    base_path = f'/Users/gwisna/Desktop/hF9mousecombined/BC{barcode}/'
    bam_path = f'{base_path}BC{barcode}_deduplicated.bam'
    samfile = pysam.AlignmentFile(bam_path, "rb")

    chromosome_name = 'ITRtoITRinsert'
    region_start = 600
    region_end = 900

    # Filter FASTQ file and get valid read names
    fastq_file = f'/Users/gwisna/Desktop/with__pear/BC{barcode}/BC{barcode}_umi.fastq'
    target_sequence = "CAAGGGAACCTTGAGAGAG"
    excluded_sequence = "GGGGGGGGGG"
    valid_read_names = filter_fastq_and_extract_read_names(fastq_file, target_sequence, excluded_sequence)
    print(f"Total valid read names: {len(valid_read_names)}") 

    # Only process reads with valid read names
    filtered_chimeric_reads, ITRtoITRinsert_reads = categorize_chimeric_reads(samfile, chromosome_name, region_start, region_end, valid_read_names)
    print(f"Total filtered chimeric reads: {len(filtered_chimeric_reads)}") 
    print(f"Total ITRtoITRinsert reads: {len(ITRtoITRinsert_reads)}")  
    # Save reads specific to 'ITRtoITRinsert'
    with open(f'{base_path}BC{barcode}_chimeric_reads_results_ITRtoITRinsert.txt', 'w') as file:
        file.write("Read Name\tOriginal Align\tChimeric Align\tCigar\n")
        for read in ITRtoITRinsert_reads:
            file.write(f"{read[0]}\t{read[1]}\t{read[2]}\t{read[3]}\n")

    samfile.close()

    # Process and analyze filtered data
    data_ITRtoITRinsert = pd.read_csv(f'{base_path}BC{barcode}_chimeric_reads_results_ITRtoITRinsert.txt', sep='\t')
    data_ITRtoITRinsert.drop_duplicates(subset='Read Name', keep='first', inplace=True)
    data_ITRtoITRinsert['Chrom'] = data_ITRtoITRinsert['Chimeric Align'].apply(lambda x: x.split(':')[0])
    data_ITRtoITRinsert['Position'] = data_ITRtoITRinsert['Chimeric Align'].apply(lambda x: int(x.split(':')[1].split(' ')[0]))
    location_counts_ITRtoITRinsert = data_ITRtoITRinsert.groupby(['Chrom', 'Position']).size().reset_index(name='Counts')
    location_counts_ITRtoITRinsert.to_csv(f'{base_path}BC{barcode}_integration_site_distribution_ITRtoITRinsert.csv', index=False)

    output_fastq = f'{base_path}BC{barcode}_chimeric_reads_results_filtered.fastq'
    filter_fastq(fastq_file, valid_read_names, output_fastq)

def filter_fastq(fastq_file, valid_read_names, output_file):
    """Filter a FASTQ file based on a set of valid read names and write the matching records to a new file."""
    with open(fastq_file, 'r') as fq, open(output_file, 'w') as out_fq:
        while True:
            header = fq.readline() 
            if not header:
                break  # End of file reached
            sequence = fq.readline().strip().upper()  
            plus = fq.readline() 
            quality = fq.readline()  

            # Extract the read name from the header 
            read_name_in_header = header.split()[0][1:] 
            if read_name_in_header in valid_read_names:
                out_fq.write(header + sequence + plus + quality)

for barcode in range(31, 41):
    process_barcode(barcode)

print("Processing complete for all files.")

