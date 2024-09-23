#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 15:55:31 2024

@author: cnelsonlab
"""

import pysam
import pandas as pd
import matplotlib.pyplot as plt

# Open the BAM file
bam_file = "/hdd/DRS/m6A/bam/IVCKM_DRS_m6a_sorted.bam"
bam = pysam.AlignmentFile(bam_file, "rb")

# Initialize lists to store data
positions = []
methylation_probs = []

# Iterate through each read in the BAM file
for read in bam.fetch():
    # Check if the read has methylation tags (MM and ML)
    if read.has_tag("MM") and read.has_tag("ML"):
        # Get methylation likelihoods
        methylation_likelihoods = read.get_tag("ML")
        
        # Get the reference position for this read
        ref_positions = read.get_reference_positions()
        
        # Save data
        for pos, prob in zip(ref_positions, methylation_likelihoods):
            positions.append(pos)
            methylation_probs.append(prob)

# Close the BAM file
bam.close()

# Convert to a DataFrame
df = pd.DataFrame({
    "Position": positions,
    "Methylation_Probability": methylation_probs
})

# Save to CSV for future reference
df.to_csv("/hdd/DRS/m6A/python_analysis/IVCKMmethylation_data.csv", index=False)

print(df.head())

##plot generation####333
df = pd.read_csv("/hdd/DRS/m6A/python_analysis/IVCKMmethylation_data.csv")

# Create a scatter plot of Methylation Probability across Genomic Position
plt.figure(figsize=(10,6))
plt.scatter(df["Position"], df["Methylation_Probability"], alpha=0.5)
plt.xlabel("Genomic Position")
plt.ylabel("Methylation Probability")
plt.title("Methylation Probability Across the Genome")
plt.show()