#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 15:30:01 2024

@author: gwisna
"""

import glob
import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import logging

# Set up logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# Input directory and file pattern
input_directory = "/Users/gwisna/Desktop/hF9mousecombined/meme/"
file_pattern = "non_aligned_sequencesBC*.fasta"
combined_output_file = "/Users/gwisna/Desktop/hF9mousecombined/meme/combined_non_aligned_sequences.fasta"
meme_output_dir = "/Users/gwisna/Desktop/hF9mousecombined/meme_output_combined"

# Function to read all FASTA files matching the pattern and combine sequences
def read_and_combine_fasta_files(directory, pattern):
    file_paths = glob.glob(directory + pattern)
    combined_sequences = []

    for file_path in file_paths:
        logging.info(f"Reading sequences from {file_path}")
        for record in SeqIO.parse(file_path, "fasta"):
            combined_sequences.append(record)

    return combined_sequences

# Step 1: Read and combine sequences from all input files
combined_sequences = read_and_combine_fasta_files(input_directory, file_pattern)

# Step 2: Write combined sequences to a new FASTA file
if combined_sequences:
    with open(combined_output_file, 'w') as output_handle:
        SeqIO.write(combined_sequences, output_handle, "fasta")
    logging.info(f"Combined sequences written to {combined_output_file}")
else:
    logging.info("No sequences found.")

# Step 3: Run MEME Suite for motif discovery on the combined 3' end sequences
def run_meme(input_fasta, output_dir, minw=6, maxw=50, nmotifs=3):
    subprocess.run([
        'meme', input_fasta, '-oc', output_dir, '-dna', '-mod', 'zoops', 
        '-nmotifs', str(nmotifs), '-minw', str(minw), '-maxw', str(maxw)
    ])

run_meme(combined_output_file, meme_output_dir, minw=6, maxw=50, nmotifs=10)

logging.info(f"MEME analysis completed. Results are saved in the {meme_output_dir} directory.")
