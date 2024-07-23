#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 08:39:01 2024

@author: harumi
"""

import subprocess
import mappy as mp
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import logging
import re

# Set up logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# Input file paths
reference_file = "/Users/gwisna/Desktop/shortF9.fa"
reads_file = "/Users/gwisna/Desktop/hF9mousecombined/BC31/BC31_chimeric_reads_results_filtered.fastq"
non_aligned_output_file = "/Users/gwisna/Desktop/hF9mousecombined/non_aligned_sequencesBC31.fasta"
meme_output_dir = "/Users/gwisna/Desktop/hF9mousecombined/meme_output"

# Function to sanitize sequence identifiers and sequences
def sanitize_identifier(identifier):
    return re.sub(r'[^a-zA-Z0-9_]', '_', identifier)

def sanitize_sequence(sequence):
    # Ensure the sequence contains only valid nucleotide characters
    return re.sub(r'[^ACGTacgt]', '', sequence)

# Function to extract non-aligning portions
def extract_non_aligning_portions(seq, alignments):
    non_aligned_regions = []
    aligned_regions = [(a.q_st, a.q_en) for a in alignments]  # Use query coordinates

    current_pos = 0
    for start, end in sorted(aligned_regions):
        if current_pos < start:
            non_aligned_regions.append((current_pos, start))
        current_pos = end
    if current_pos < len(seq):
        non_aligned_regions.append((current_pos, len(seq)))

    non_aligned_sequences = [seq[start:end] for start, end in non_aligned_regions]
    return non_aligned_sequences

# Step 1: Load the reference sequence
logging.info(f"Loading reference sequence from {reference_file}")
ref_seq = mp.Aligner(reference_file, preset='sr')  # Using 'sr' preset for short reads

if not ref_seq:
    raise Exception("Failed to load/build index")

# Step 2: Align reads and collect non-aligning portions
non_aligned_sequences = []

try:
    logging.info(f"Reading FASTQ file from {reads_file}")
    for name, seq, qual in mp.fastx_read(reads_file):
        logging.debug(f"Read name: {name}")
        logging.debug(f"Sequence: {seq}")
        alignments = list(ref_seq.map(seq))
        
        if alignments:
            logging.debug(f"Alignments found: {alignments}")
            non_aligned_seqs = extract_non_aligning_portions(seq, alignments)
            for non_aligned_seq in non_aligned_seqs:
                sanitized_id = sanitize_identifier(name)
                sanitized_seq = sanitize_sequence(non_aligned_seq)
                non_aligned_sequences.append(SeqRecord(Seq(sanitized_seq), id=sanitized_id, description="non-aligned"))
        else:
            sanitized_id = sanitize_identifier(name)
            sanitized_seq = sanitize_sequence(seq)
            non_aligned_sequences.append(SeqRecord(Seq(sanitized_seq), id=sanitized_id, description="completely non-aligned"))
except ValueError as e:
    logging.error(f"Error reading FASTQ file: {e}")

# Step 3: Write non-aligning sequences to a FASTA file
if non_aligned_sequences:
    with open(non_aligned_output_file, 'w') as output_handle:
        SeqIO.write(non_aligned_sequences, output_handle, "fasta")
    logging.info(f"Non-aligned sequences written to {non_aligned_output_file}")
else:
    logging.info("No non-aligned sequences found.")
    
# Function to extract the last 30-50 base pairs of a sequence
def extract_3_end(seq, length=50):
    return seq[-length:]

# Step 4: Run MEME Suite for motif discovery
def run_meme(input_fasta, output_dir):
    # Read the sequences from the input FASTA file
    records = list(SeqIO.parse(input_fasta, "fasta"))
    
    # Extract the 3' end of each sequence (30-50 bp)
    end_sequences = []
    for record in records:
        if len(record.seq) >= 30:  # Ensure the sequence is at least 30 bp long
            end_seq = record.seq[-50:]  # Extract up to the last 50 bp
            end_sequences.append(SeqRecord(end_seq, id=record.id, description="3' end non-aligned"))

    # Write the 3' end sequences to a new FASTA file
    end_fasta = input_fasta.replace(".fasta", "_3end.fasta")
    with open(end_fasta, 'w') as output_handle:
        SeqIO.write(end_sequences, output_handle, "fasta")
    
    # Run MEME on the 3' end sequences
    subprocess.run(['meme', end_fasta, '-oc', output_dir, '-dna', '-mod', 'zoops', '-nmotifs', '10', '-minw', '6', '-maxw', '50'])

run_meme(non_aligned_output_file, meme_output_dir)

# Step 5: Run MEME Suite for motif discovery on the last 30-50 bp of the sequences
def run_meme(input_fasta, output_dir, minw=6, maxw=50, nmotifs=3):
    # Read the sequences from the input FASTA file
    records = list(SeqIO.parse(input_fasta, "fasta"))
    
    # Extract the 3' end of each sequence
    end_sequences = []
    for record in records:
        if len(record.seq) >= 30:  # Ensure the sequence is at least 30 bp long
            end_seq = extract_3_end(record.seq, 50)  # Extract up to 50 bp
            end_sequences.append(SeqRecord(end_seq, id=record.id, description="3' end non-aligned"))

    # Write the 3' end sequences to a new FASTA file
    end_fasta = input_fasta.replace(".fasta", "_3end.fasta")
    with open(end_fasta, 'w') as output_handle:
        SeqIO.write(end_sequences, output_handle, "fasta")
    
    # Run MEME on the 3' end sequences with adjusted parameters
    subprocess.run([
        'meme', end_fasta, '-oc', output_dir, '-dna', '-mod', 'zoops', 
        '-nmotifs', str(nmotifs), '-minw', str(minw), '-maxw', str(maxw)
    ])

run_meme(non_aligned_output_file, meme_output_dir, minw=6, maxw=50, nmotifs=5)
logging.info(f"MEME analysis completed. Results are saved in the {meme_output_dir} directory.")