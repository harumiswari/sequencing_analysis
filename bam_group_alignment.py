#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 14:55:06 2024

@author: gwisna
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 11:12:56 2024

@author: gwisna
"""

import pysam
import mappy as mp
import matplotlib.pyplot as plt
import os
from Bio.Seq import Seq
import pandas as pd

input_dir = "/Users/gwisna/Desktop/Ckm_mRNA_picard/"
base_file_name = "BC"
barcode_numbers = list(range(21, 30))

# Dictionary to store index paths and their respective expected sequences
indexes = {
    "/Users/gwisna/Desktop/Ckm_mRNA.fa": "ACCTCCACAGCACAGACAGACACTCAGGAGCCAGCCAGCCAG",
    "/Users/gwisna/Desktop/F9mRNARACE.fa": "ACCTCCACAGCACAGACAGACACTCAGGAG"
}
#indexes ={"/Users/gwisna/Desktop/tn5_DNAF9.fa":"TGTGACTGTCCTTTTCTGGCTTTAGACAAAAGGTTTTGCC" }

# Parameters
min_read_length = 75
results = []
coverage_results = []

# Process each BAM file for each index and expected sequence
for index_path, expected_sequence in indexes.items():
    expected_sequence_rc = str(Seq(expected_sequence).reverse_complement())

    # Load the index
    index = mp.Aligner(index_path)
    if not index:
        raise Exception(f"Index file {index_path} not found")

    for barcode in barcode_numbers:
        bam_file = os.path.join(input_dir, f"{base_file_name}{barcode}/{base_file_name}{barcode}.RG_rmdup_sort.bam")
        
        # Check if the input BAM file exists
        if not os.path.exists(bam_file):
            print(f"Input file {bam_file} not found, skipping")
            continue

        # Initialize coverage alignment list
        #coverage_alignment = [0] * 1451
        coverage_alignment = [0] * 3386

        # Total number of mapped reads
        total_mapped_reads = 0

        # Open the BAM file
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            # Count total number of mapped reads
            for read in bam:
                if not read.is_unmapped and \
                    len(read.query_sequence) > min_read_length and \
                    (expected_sequence in read.query_sequence or expected_sequence_rc in read.query_sequence):
                        total_mapped_reads += 1

        print(f"Total number of mapped reads in {bam_file} with index {os.path.basename(index_path)}: {total_mapped_reads}")

        with pysam.AlignmentFile(bam_file, "rb") as bam:
            # Process each read in the BAM file
            for read in bam:
                if not read.is_unmapped and \
                    len(read.query_sequence) > min_read_length and \
                    (expected_sequence in read.query_sequence or expected_sequence_rc in read.query_sequence):
                        for hit in index.map(read.seq):  
                            if hit.mapq > 30:
                                for i in range(hit.r_st, hit.r_en):
                                    coverage_alignment[i] += 1
        
        # Coverage results
        for barcode in barcode_numbers:
            for position, coverage in enumerate(coverage_alignment):
                coverage_results.append({
                    'Barcode': f"{base_file_name}{barcode}",
                    'Position': position,
                    'Coverage': coverage
                })

        # Plot coverage alignment
        plt.figure(figsize=(10, 6))
        plt.plot(coverage_alignment)
        plt.xlabel('Position')
        plt.ylabel('Coverage')
        plt.title(f'Coverage Alignment for {base_file_name}{barcode} with {os.path.basename(index_path)}')
        plt.show()

        alignment_counts = {}

        # Process each read in the BAM file
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam:
                if not read.is_unmapped and len(read.query_sequence) > min_read_length and \
                    (expected_sequence in read.query_sequence or expected_sequence_rc in read.query_sequence):
                    for hit in index.map(read.seq):  
                        if hit.mapq > 30:
                            ref_name = hit.ctg
                            if ref_name not in alignment_counts:
                                alignment_counts[ref_name] = {'forward': 0, 'reverse': 0}
                            if read.is_reverse:
                                alignment_counts[ref_name]['reverse'] += 1
                            else:
                                alignment_counts[ref_name]['forward'] += 1

        # Alignment count results
        for ref_name, counts in alignment_counts.items():
            results.append({
                'Barcode': f"{base_file_name}{barcode}",
                'Reference sequence': ref_name,
                'Forward strand': counts['forward'],
                'Reverse strand': counts['reverse'],
                'Total mapped reads': total_mapped_reads
            })


df = pd.DataFrame(results)
df_coverage = pd.DataFrame(coverage_results)


print(df)
print(df_coverage)

input_dir_name = os.path.basename(os.path.normpath(input_dir))
csv_file_path = os.path.join(input_dir, f"{input_dir_name}_BC1-20_alignment_counts.csv")
df.to_csv(csv_file_path, index=False)

# To save the alignment coverage in each barcode
for barcode in barcode_numbers:
    # Filter DataFrame for the current barcode
    df_barcode = df_coverage[df_coverage['Barcode'] == f"{base_file_name}{barcode}"]

    # Save DataFrame to a CSV file for the current barcode
    csv_file_path = f"{input_dir}{base_file_name}{barcode}.csv"
    df_barcode.to_csv(csv_file_path, index=False)
