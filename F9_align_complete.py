#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 11:12:56 2024

@author: harumi
"""

import pysam
import mappy as mp
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import os
import pandas as pd

input_dir = "/Users/gwisna/Desktop/with__pear/"
index_path = "/Users/gwisna/Desktop/F9mRNARACE.fa"
base_file_name = "BC"
barcode_numbers = list(range(21, 31))

# Parameters
min_read_length = 125
expected_sequence = "ACCTCCACAGCACAGACAGACACTCAGGAG"
expected_sequence_rc = str(Seq(expected_sequence).reverse_complement())
max_consecutive_gs = 6  # Threshold for consecutive 'G's considered as low quality

bam_file_template = os.path.join(input_dir, f"{base_file_name}{{barcode}}/{base_file_name}{{barcode}}_deduplicated.bam")

def is_low_quality(read_sequence, max_consecutive_gs):
    """
    Check if the read sequence has a stretch of consecutive Gs exceeding the threshold.
    """
    return 'G' * max_consecutive_gs in read_sequence

results = []
coverage_results = []

# Process that each Bam file will go through (alignment, fastq write, counting)
def process_bam_file(bam_file, index_path, barcode):
    total_mapped_reads = 0
    coverage_alignment = [0] * 1266 #adjust this
    alignment_counts = {'forward': 0, 'reverse': 0}

    index = mp.Aligner(index_path)
    if not index:
        raise Exception(f"Index file {index_path} not found")

    # Open the BAM file and the output fastq file
    output_fastq_file = os.path.join(input_dir, f"F9_aligned_reads_BC{barcode}.fastq")
    with pysam.AlignmentFile(bam_file, "rb") as bam, open(output_fastq_file, "w") as fastq_out:
        # Count total number of mapped reads
        for read in bam:
            if not read.is_unmapped and \
                len(read.query_sequence) > min_read_length and \
                (expected_sequence in read.query_sequence or expected_sequence_rc in read.query_sequence) and \
                not is_low_quality(read.query_sequence, max_consecutive_gs):
                
                for hit in index.map(read.seq):  
                    if hit.mapq > 30:
                        total_mapped_reads += 1
                        #write aligned reads to fastq
                        fastq_out.write(f"@{read.query_name}\n")
                        fastq_out.write(f"{read.query_sequence}\n")
                        fastq_out.write("+\n")
                        fastq_out.write(f"{read.qual}\n")
                        
                        # Count coverage alignment based on reference position
                        for i in range(hit.r_st, hit.r_en):
                            if i < 1266:  #adjust this
                                coverage_alignment[i] += 1
                                
                        # Count alignment by strand orientation
                        if read.is_reverse:
                            alignment_counts['reverse'] += 1
                        else:
                            alignment_counts['forward'] += 1

    # Plot coverage alignment
    plt.figure(figsize=(10, 6))
    plt.plot(coverage_alignment)
    plt.xlabel('Position')
    plt.ylabel('Coverage')
    plt.title(f'Coverage Alignment for {base_file_name}{barcode}')
    plt.show()

    # Append coverage alignment results for this barcode
    for position, coverage in enumerate(coverage_alignment):
        coverage_results.append({
            'Barcode': f"{base_file_name}{barcode}",
            'Position': position,
            'Coverage': coverage
        })

    # Append alignment count results for this barcode
    for ref_name, counts in alignment_counts.items():
        results.append({
            'Barcode': f"{base_file_name}{barcode}",
            'Forward strand': alignment_counts['forward'],
            'Reverse strand': alignment_counts['reverse'],
            'Total mapped reads': total_mapped_reads
        })

    return total_mapped_reads

# Process each barcode
for barcode in barcode_numbers:
    bam_file = bam_file_template.format(barcode=barcode)

    if not os.path.exists(bam_file):
        print(f"BAM file {bam_file} not found, skipping barcode {barcode}")
        continue

    total_mapped_reads = process_bam_file(bam_file, index_path, barcode)
    print(f"Processed {base_file_name}{barcode}: Total mapped reads = {total_mapped_reads}")

df = pd.DataFrame(results)
df_coverage = pd.DataFrame(coverage_results)

#CSV for summary of alignment
csv_file_path = os.path.join(input_dir, f"{os.path.basename(os.path.normpath(input_dir))}F9_alignment_counts.csv")
df.to_csv(csv_file_path, index=False)

#CSV for individual barcode
for barcode in barcode_numbers:
    df_barcode_coverage = df_coverage[df_coverage['Barcode'] == f"{base_file_name}{barcode}"]
    csv_file_path_coverage = os.path.join(input_dir, f"{base_file_name}{barcode}F9_coverage_counts.csv")
    df_barcode_coverage.to_csv(csv_file_path_coverage, index=False)

print("Alignment Counts:")
print(df)
print("\nCoverage Counts:")
print(df_coverage)
