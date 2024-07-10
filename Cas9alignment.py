#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 30 16:11:30 2024

@author: gwisna
"""

import pandas as pd
import pysam
import mappy as mp


def categorize_all_reads(samfile, target_chrom):
    all_reads = []
    read_names_processed = set()

    for read in samfile.fetch(target_chrom):
        if read.query_name in read_names_processed or read.is_unmapped:
            continue
        read_names_processed.add(read.query_name)
        original_chrom = read.reference_name
        original_pos = read.reference_start
        end_pos = original_pos + read.reference_length  # Calculate the end position of the read alignment
        cigar = read.cigarstring

        # Only include reads that align between 622 and 4857 bp or those that align within these ranges
        if (original_pos >= 622 and original_pos <= 4857) or (end_pos >= 622 and end_pos <= 4857):
            read_details = (read.query_name, f"{original_chrom}:{original_pos}", cigar)
            all_reads.append(read_details)

    return all_reads

def process_barcode(barcode):
    base_path = f'/Users/gwisna/Desktop/ITRCas9_nopear/BC{barcode}/'
    bam_path = f'{base_path}BC{barcode}_deduplicated.bam'
    samfile = pysam.AlignmentFile(bam_path, "rb")

    target_chrom = 'ITRtoITRCas9'

    all_reads = categorize_all_reads(samfile, target_chrom)

    # Convert to DataFrame and deduplicate
    all_reads_df = pd.DataFrame(all_reads, columns=['Read Name', 'Original Align', 'Cigar'])
    all_reads_df.drop_duplicates(subset='Read Name', keep='first', inplace=True)

    # Extract chromosome and position
    all_reads_df['Chrom'] = all_reads_df['Original Align'].apply(lambda x: x.split(':')[0])
    all_reads_df['Position'] = all_reads_df['Original Align'].apply(lambda x: int(x.split(':')[1].split(' ')[0]))

    # Process data
    data = all_reads_df[all_reads_df['Chrom'] == target_chrom]
    location_counts = data.groupby(['Chrom', 'Position']).size().reset_index(name='Counts')
    output_path_counts = f'{base_path}BC{barcode}_{target_chrom}_reads_counts.csv'
    location_counts.to_csv(output_path_counts, index=False)
    data[['Read Name']].to_csv(f'{base_path}BC{barcode}_{target_chrom}_read_names.txt', index=False, header=False)

    samfile.close()

    print(f"Read name lists and counts from {target_chrom} saved to {output_path_counts}")

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
    data_path = f'/Users/gwisna/Desktop/ITRCas9_nopear/BC{barcode}/BC{barcode}_ITRtoITRCas9_read_names.txt'
    data = pd.read_csv(data_path, delimiter='\t', header=None, names=['Read Name'])
    read_names = set(data['Read Name'])

    fastq_file = f'/Users/gwisna/Desktop/with__pear/BC{barcode}/BC{barcode}_umi.fastq'
    output_file = f'/Users/gwisna/Desktop/ITRCas9_nopear/BC{barcode}/BC{barcode}_ITRtoITRCas9.fastq'

    filter_fastq(fastq_file, read_names, output_file)


def align_and_filter_reads_mappy(barcode):
    base_path = f'/Users/gwisna/Desktop/ITRCas9_nopear/BC{barcode}/'
    fastq_file = f'{base_path}BC{barcode}_ITRtoITRCas9.fastq'
    reference_fasta = '/Users/gwisna/Desktop/ITRtoITRCas9.fa'
    filtered_fastq = f'{base_path}BC{barcode}_filtered.fastq'
    output_sam = f'{base_path}BC{barcode}_filtered.sam'
    output_bam = f'{base_path}BC{barcode}_filtered.bam'
    sorted_bam = f'{base_path}BC{barcode}_filtered_sorted.bam'
    coverage_output = f'{base_path}BC{barcode}_coverage.txt'

    # Initialize the aligner
    aligner = mp.Aligner(reference_fasta)
    if not aligner:
        raise Exception("ERROR: failed to load/build index")

    read_count = 0
    has_alignment = False

    # Read the FASTQ file and align
    with open(fastq_file, 'r') as fq, open(filtered_fastq, 'w') as out_fq, open(output_sam, 'w') as out_sam:
        while True:
            header = fq.readline()
            if not header:
                break
            sequence = fq.readline().strip()
            plus = fq.readline()
            quality = fq.readline().strip()

            # Align the read
            alignments = aligner.map(sequence)

            # Check if the read matches 100% with the reference
            match_100_percent = False
            for aln in alignments:
                if aln.mlen == len(sequence) and aln.NM == 0:  # 100% match
                    match_100_percent = True
                    break

            # Keep the read if it does not match 100%
            if not match_100_percent:
                out_fq.write(f"{header}{sequence}\n{plus}{quality}\n")
                read_count += 1
                has_alignment = True
                # Write SAM alignment
                for aln in alignments:
                    if aln.is_primary:
                        out_sam.write(aln.to_sam(header, sequence, quality))
                        break

    # Save the read count to a file
    read_count_path = f'{base_path}BC{barcode}_filtered_read_count.txt'
    with open(read_count_path, 'w') as count_file:
        count_file.write(f'Read Count: {read_count}\n')

    print(f"Filtered reads saved to {filtered_fastq} with {read_count} reads")

    if not has_alignment:
        print(f"No valid alignments found for {barcode}. Skipping BAM and coverage generation.")
        return

    # Convert SAM to BAM, sort and index using pysam
    with pysam.AlignmentFile(output_sam, "r") as samfile:
        with pysam.AlignmentFile(output_bam, "wb", header=samfile.header) as bamfile:
            for s in samfile:
                bamfile.write(s)
    
    pysam.sort("-o", sorted_bam, output_bam)
    pysam.index(sorted_bam)

    # Generate coverage
    coverage = pysam.depth("-aa", sorted_bam)
    with open(coverage_output, 'w') as cov_out:
        cov_out.write(coverage)

    print(f"Coverage map saved to {coverage_output}")

# Loop over barcodes
for barcode in range(11, 21):
    align_and_filter_reads_mappy(barcode)