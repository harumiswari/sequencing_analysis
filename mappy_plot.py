#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 23 10:16:22 2024

@author: gwisna
"""

import mappy as mp
import matplotlib.pyplot as plt

# Load the reference genome
reference = "/Users/gwisna/Desktop/ITRtoITRinsert.fa"
aligner = mp.Aligner(reference)
if not aligner:
    raise ValueError("ERROR: failed to load/build index")

# Define regions in the fasta file
regions = {
    'ITR1': (0, 146),
    'U6_gRNA': (147, 387),
    'F9': (388, 3900),
    'ITR2': (3901, 4583)
}

# Define the barcodes and corresponding colors
barcodes = [5,10,15,20]
colors = ['blue', 'green', 'red', 'purple']

# Dictionary to store counts of forward and reverse reads for each region for each barcode
barcode_counts = {barcode: {region: {'forward': 0, 'reverse': 0} for region in regions} for barcode in barcodes}

def align_and_count(fastq_file, barcode, mapq_threshold=20):
    read_processed = set()
    for name, seq, qual in mp.fastx_read(fastq_file):
        if name in read_processed:
            continue
        has_alignment = False
        for hit in aligner.map(seq):
            if hit.mapq >= mapq_threshold and not has_alignment:
                # Determine the region
                for region_name, (start, end) in regions.items():
                    if hit.r_st >= start and hit.r_en <= end:
                        direction = 'forward' if hit.strand == 1 else 'reverse'
                        barcode_counts[barcode][region_name][direction] += 1
                        break
                has_alignment = True
        read_processed.add(name)

# Process each FASTQ file corresponding to each barcode
for barcode in barcodes:
    fastq_path = f"/Users/gwisna/Desktop/hF9mousecombined/BC{barcode}/BC{barcode}_ITRtoITRinsert.fastq"
    align_and_count(fastq_path, barcode)

# Output results
for barcode in barcodes:
    with open(f'/Users/gwisna/Desktop/hF9mousecombined/BC{barcode}/BC{barcode}_region_orientation_counts.txt', 'w') as f:
        f.write(f"Orientation and Region Results for Barcode {barcode}:\n")
        for region in regions:
            f.write(f"Region {region}:\n")
            for direction in ['forward', 'reverse']:
                count = barcode_counts[barcode][region][direction]
                f.write(f"  {direction.title()} Reads: {count}\n")
