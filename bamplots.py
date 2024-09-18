#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 20:17:31 2024

@author: gwisna
"""

import pysam
import matplotlib.pyplot as plt
import numpy as np

# Path to your BAM file and synthetic reference FASTA file
bam_file = "/hdd/Long_read_alignment/Ckmgene_alignment/directdnaivck2secondrun_Ckmmouse_sorted.bam"
synthetic_fasta = "/hdd/DRS/Ref_file/Ckmmouse.fa"

# Open the BAM file with pysam
samfile = pysam.AlignmentFile(bam_file, "rb")

# Load the synthetic reference FASTA file
ref_fasta = pysam.FastaFile(synthetic_fasta)

# Specify the synthetic reference sequence name and region of interest
synthetic_sequence_name = "NC_000073.7:19145019-19155508"  # Example name in your synthetic FASTA
start_position = 1  # Change to your desired start position
end_position = 10490   # Change to your desired end position

# Initialize an array for coverage data
coverage = np.zeros(end_position - start_position)

# Generate pileup for the specified region
for pileupcolumn in samfile.pileup(synthetic_sequence_name, start_position, end_position):
    pos = pileupcolumn.reference_pos - start_position
    if 0 <= pos < len(coverage):
        coverage[pos] = pileupcolumn.nsegments  # Number of reads covering this position

# Close the BAM file
samfile.close()

# Optionally, fetch the reference nucleotide sequence from the synthetic FASTA
reference_sequence_nucleotides = ref_fasta.fetch(synthetic_sequence_name, start_position, end_position)
print(reference_sequence_nucleotides)  # Display the synthetic reference sequence in the specified range

# Close the FASTA file
ref_fasta.close()

# Plot the coverage
plt.figure(figsize=(10, 4))
plt.plot(range(start_position, end_position), coverage, color='blue')
plt.title(f"Read Pileup for {synthetic_sequence_name}:{start_position}-{end_position}")
plt.xlabel("Position")
plt.ylabel("Coverage")
plt.grid(True)
#plt.savefig('/Users/gwisna/Library/CloudStorage/Box-Box/NELSONLAB/Direct RNA Seq/IVT_TrialRun_1/EGFP/IVT_GFP_read_pileup.svg', format='svg')
plt.show()

###############################################3

import pysam
import matplotlib.pyplot as plt

# Path to your BAM file and synthetic reference FASTA file
bam_file = "/hdd/LR_assembly_flye/directdnaseqIVCK2_medaka_output/calls_to_draft.bam"
samfile = pysam.AlignmentFile(bam_file, "rb")
print(samfile.header)

synthetic_sequence_name = "NC_000073.7:19145019-19155508"
start_position = 1  # Change to your desired start position
#end_position = 4583  # Change to your desired end position
end_position = 80000

# Open the BAM file with pysam
samfile = pysam.AlignmentFile(bam_file, "rb")

# Initialize lists to hold the start and end positions of reads
read_starts = []
read_ends = []

# Initialize a variable to keep track of the y-position for each read
y_positions = []

# Iterate over each read in the BAM file for the specified region
for i, read in enumerate(samfile.fetch(synthetic_sequence_name, start_position, end_position)):
    read_starts.append(read.reference_start)
    read_ends.append(read.reference_end)
    y_positions.append(i)  # Use the index `i` to set a unique y-position for each read

# Close the BAM file
samfile.close()

# Create the figure for the plot
plt.figure(figsize=(10, 6))

# Plot each read as a line from its start to end position with a unique y-position
for start, end, y in zip(read_starts, read_ends, y_positions):
    plt.plot([start, end], [y, y], color='blue', alpha=0.5)
# **Zoom in on the region of interest (2 kb - 6 kb)**
plt.xlim(0, 10000)  # Adjust the x-axis limits to zoom in on 2kb-6kb
plt.ylim(-1, 1500)  # Adjust y-limits to fit all reads
plt.title(f"Individual Read Coverage for {synthetic_sequence_name}:2000-6000")
plt.xlabel("Position")
plt.ylabel("Reads")
plt.grid(True)
plt.show()


# Adjust plot limits and labels
plt.xlim(start_position, end_position)
plt.ylim(-1, len(y_positions))  # Adjust y-limits to fit all reads
plt.title(f"Individual Read Coverage for{synthetic_sequence_name}:{start_position}-{end_position}")
plt.xlabel("Position")
plt.ylabel("Reads")
plt.grid(True)
#plt.savefig('/Users/gwisna/Library/CloudStorage/Box-Box/NELSONLAB/Direct RNA Seq/IVT_TrialRun_1/EGFP/IVT_GFP_individual_read_coverage.svg', format='svg')
plt.show()


############################################

import pysam
import matplotlib.pyplot as plt

# Path to your BAM file and synthetic reference FASTA file
bam_file = "/Users/gwisna/Library/CloudStorage/Box-Box/NELSONLAB/Direct RNA Seq/DRS-LentiGFP/DRS_lenti_GFP_sorted.bam"
synthetic_sequence_name = "Untitled"  # Example name in your synthetic FASTA
start_position = 1  # Change to your desired start position
end_position = 5337  # Change to your desired end position

# Open the BAM file with pysam
samfile = pysam.AlignmentFile(bam_file, "rb")

# Initialize lists to hold the start and end positions of reads
read_starts = []
read_ends = []

# Initialize a variable to keep track of the y-position for each read
y_positions = []

# Iterate over each read in the BAM file for the specified region
for i, read in enumerate(samfile.fetch(synthetic_sequence_name, start_position, end_position)):
    read_starts.append(read.reference_start)
    read_ends.append(read.reference_end)
    y_positions.append(i)  # Use the index `i` to set a unique y-position for each read

# Close the BAM file
samfile.close()

# Create the figure for the plot
plt.figure(figsize=(10, 6))

# Plot each read as a line from its start to end position with a unique y-position
for start, end, y in zip(read_starts, read_ends, y_positions):
    plt.plot([start, end], [y, y], color='blue', alpha=0.5)

# Adjust plot limits and labels
plt.xlim(start_position, end_position)
plt.ylim(-1, len(y_positions))  # Adjust y-limits to fit all reads
plt.title(f"Individual Read Coverage for Lenti-GFP:{start_position}-{end_position}")
plt.xlabel("Position")
plt.ylabel("Reads")
plt.grid(True)
plt.savefig('/Users/gwisna/Library/CloudStorage/Box-Box/NELSONLAB/Direct RNA Seq/DRS-LentiGFP/DRS_lenti_GFP_individual_read_coverage.svg', format='svg')
plt.show()

######################################

import pysam

# Path to your BAM file
bam_file = "/Users/gwisna/Library/CloudStorage/Box-Box/NELSONLAB/Direct RNA Seq/DRS-LentiGFP/DRS_lenti_GFP_sorted.bam"

# Open the BAM file with pysam
samfile = pysam.AlignmentFile(bam_file, "rb")

# Initialize a list to store read lengths and read information
read_lengths = []

# Iterate over each read in the BAM file
for read in samfile.fetch():
    read_length = read.query_length  # Get the length of the read
    read_lengths.append((read.query_name, read_length))  # Store the read name and length

# Close the BAM file
samfile.close()

# Sort the reads by length in descending order
read_lengths.sort(key=lambda x: x[1], reverse=True)

# Print the top N longest reads (e.g., top 10 longest reads)
top_n = 10
print(f"Top {top_n} longest reads:")
for read_name, read_length in read_lengths[:top_n]:
    print(f"Read Name: {read_name}, Length: {read_length}")
    
###################################################3

import pysam
import matplotlib.pyplot as plt

# Path to your BAM file
bam_file = "/Users/gwisna/Library/CloudStorage/Box-Box/NELSONLAB/Direct RNA Seq/DRS-LentiGFP/DRS_lenti_GFP_sorted.bam"

# Open the BAM file with pysam
samfile = pysam.AlignmentFile(bam_file, "rb")

# Initialize a list to store mapping quality scores
mapping_qualities = []

# Iterate over each read in the BAM file
for read in samfile.fetch():
    if not read.is_unmapped:  # Only consider mapped reads
        mapping_qualities.append(read.mapping_quality)  # Extract the mapping quality

# Close the BAM file
samfile.close()

# Plot the distribution of mapping qualities
plt.figure(figsize=(10, 6))
plt.hist(mapping_qualities, bins=50, color='blue', alpha=0.7, edgecolor='black')
plt.title('Alignment Quality Distribution')
plt.xlabel('Mapping Quality Score')
plt.ylabel('Frequency')
plt.grid(False)
plt.savefig('/Users/gwisna/Library/CloudStorage/Box-Box/NELSONLAB/Direct RNA Seq/DRS-LentiGFP/DRS_lenti_GFP_alignment_quality_distribution.svg', format='svg')
plt.show()

###############################################

import pysam
import matplotlib.pyplot as plt

# Path to your BAM file and consensus FASTA file
bam_file = "/hdd/LR_assembly_flye/directdnaseqIVCK2_medaka_output/calls_to_draft.bam"
samfile = pysam.AlignmentFile(bam_file, "rb")

# Print BAM header to inspect available references (contigs)
print(samfile.header)

# Set start and end positions (you can adjust this as necessary)
start_position = 1
end_position = 80000  # Adjust based on your region of interest

# Initialize a dictionary to store read positions for each contig
contig_read_positions = {}

# Iterate over each contig (reference sequence) in the BAM file
for contig in samfile.references:
    # Initialize lists to hold the start and end positions of reads for the current contig
    read_starts = []
    read_ends = []
    
    # Iterate over each read in the BAM file for the specified contig
    for i, read in enumerate(samfile.fetch(contig, start_position, end_position)):
        read_starts.append(read.reference_start)
        read_ends.append(read.reference_end)
    
    # Store the start and end positions for the contig
    contig_read_positions[contig] = (read_starts, read_ends)

# Close the BAM file
samfile.close()

# Plot the read coverage for each contig
plt.figure(figsize=(10, 6))

# Iterate over the contig_read_positions dictionary to plot each contig's reads
for contig, (read_starts, read_ends) in contig_read_positions.items():
    if read_starts:  # Plot only if there are reads in this contig
        y_positions = range(len(read_starts))
        plt.plot([read_starts, read_ends], [y_positions, y_positions], color="blue", label=f"Contig: {contig}")

# Customize the plot
plt.xlabel("Position")
plt.ylabel("Read Index")
plt.title("BAM Read Coverage across Contigs")
plt.legend()
plt.savefig('/hdd/LR_assembly_flye/medaka.png', format='png')

