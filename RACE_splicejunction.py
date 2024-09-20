#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 15:44:05 2024

@author: cnelsonlab
"""
#############for multi track view#############

import pysam
import matplotlib.pyplot as plt
import numpy as np

def parse_gtf(gtf_file, region):
    exons = []
    chrom, region_range = region.split(":")
    region_start, region_end = map(int, region_range.split("-"))
    
    gtf = pd.read_csv(gtf_file, sep="\t", comment="#", header=None, names=[
        "seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
    
    for _, row in gtf.iterrows():
        if row['seqname'] == chrom and row['feature'] == "exon":
            if (row['start'] >= region_start and row['end'] <= region_end):
                exons.append((row['start'], row['end']))
    
    return exons

def plot_combined_coverage_and_junctions(bam_files, gtf_file, region, output_file):
    """
    Generate a combined plot of coverage across multiple samples.
    """
    fig, axes = plt.subplots(len(bam_files), figsize=(10, 6*len(bam_files)), sharex=True)
    
    chrom, region_range = region.split(":")
    region_start, region_end = map(int, region_range.split("-"))
    
    exon_colors = {
        (1, 53): 'green',        # Exon 1 (Ckm)
        (54, 59): 'blue',        # Exon 2 (Intermediary)
        (60, 1266): 'purple'     # Exon 3 (F9)
    }

    exons = parse_gtf(gtf_file, region)

    # Iterate over BAM files and plot each in a separate track
    for i, bam_file in enumerate(bam_files):
        ax = axes[i]
        junctions = {}
        coverage = np.zeros(region_end - region_start + 1)

        # Process BAM file
        samfile = pysam.AlignmentFile(bam_file, "rb")
        for read in samfile.fetch(chrom, region_start, region_end):
            if read.is_unmapped or read.is_secondary:
                continue

            # Update coverage
            for pos in range(max(region_start, read.reference_start), min(region_end, read.reference_end)):
                coverage[pos - region_start] += 1

            # Process splice junctions (CIGAR N operation)
            for cigar in read.cigartuples:
                if cigar[0] == 3:  # CIGAR N operation: skip region (splice junction)
                    start = read.reference_start
                    end = start + cigar[1]
                    if (start, end) not in junctions:
                        junctions[(start, end)] = 0
                    junctions[(start, end)] += 1

        # Plot read coverage
        ax.fill_between(np.arange(region_start, region_end + 1), coverage, step="mid", color="orange", alpha=0.5)

        # Plot junctions (arcs)
        for (start, end), count in junctions.items():
            arc_x = np.linspace(start, end, 100)
            arc_y = count - (count * (arc_x - start) * (arc_x - end) / ((end - start)**2))  # Parabolic arc
            ax.plot(arc_x, arc_y, color="red", lw=count / 5.0)  # Thickness based on read count
            ax.text((start + end) / 2, count + 2, f"{count}", fontsize=10, color="black", ha="center")  # Label with read count

        # Plot exons at the bottom
        for exon_start, exon_end in exons:
            color = exon_colors.get((exon_start, exon_end), 'gray')
            ax.add_patch(plt.Rectangle((exon_start, 0), exon_end - exon_start, 1, color=color, alpha=0.7))

        # Set labels and titles for each track
        ax.set_ylabel(f"Sample {i+1}")
        ax.set_xlim([region_start, region_end])
        ax.set_ylim([0, max(coverage) + 20])  # Adjust y-limits for better view of arcs
        ax.set_title(f"Track {i+1}: {bam_file.split('/')[-1]}")

    # Set shared x-axis label for genomic position
    axes[-1].set_xlabel(f"Genomic Position: {region}")
    
    # Save plot as SVG
    plt.savefig(output_file, format="svg")
    plt.show()

#bam_files = [
#    "/hdd/InvivoRACE__concatenated_fastq/Invivo_RACE_alignment_redo/barcode08_sorted.bam",
#    "/hdd/InvivoRACE__concatenated_fastq/Invivo_RACE_alignment_redo/barcode09_sorted.bam",
#    "/hdd/InvivoRACE__concatenated_fastq/Invivo_RACE_alignment_redo/barcode10_sorted.bam",
#    "/hdd/InvivoRACE__concatenated_fastq/Invivo_RACE_alignment_redo/barcode11_sorted.bam"
#]

bam_files = [
    "/hdd/STAR/Output/Sample_5_Aligned.sortedByCoord.out.bam",
    "/hdd/STAR/Output/Sample_6_Aligned.sortedByCoord.out.bam",
    "/hdd/STAR/Output/Sample_7_Aligned.sortedByCoord.out.bam",
    "/hdd/STAR/Output/Sample_8_Aligned.sortedByCoord.out.bam",
    "/hdd/STAR/Output/Sample_13_Aligned.sortedByCoord.out.bam",
    "/hdd/STAR/Output/Sample_14_Aligned.sortedByCoord.out.bam",
    "/hdd/STAR/Output/Sample_15_Aligned.sortedByCoord.out.bam",
#    "/hdd/STAR/Output/Sample_12_Aligned.sortedByCoord.out.bam",
]


gtf_file = "/home/cnelsonlab/F9mRNARACE.gtf"
region = "F9mRNARACE:1-1266"  # Define the region of interest
output_file = "scrambled_shortread_multiple_tract_coverage.svg"

plot_combined_coverage_and_junctions(bam_files, gtf_file, region, output_file)


################RACE_splice junction plot###################
import pysam
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors

def parse_gtf(gtf_file, region):
    """
    Parse the GTF file to extract exon information for the region of interest.
    """
    exons = []
    chrom, region_range = region.split(":")
    region_start, region_end = map(int, region_range.split("-"))
    
    gtf = pd.read_csv(gtf_file, sep="\t", comment="#", header=None, names=[
        "seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
    
    for _, row in gtf.iterrows():
        if row['seqname'] == chrom and row['feature'] == "exon":
            if (row['start'] >= region_start and row['end'] <= region_end):
                exons.append((row['start'], row['end']))
    
    return exons

def plot_sashimi_with_actual_read_counts(bam_files, gtf_file, region, output_file):
    """
    Generate a Sashimi plot from multiple BAM files and a GTF annotation file.
    Arcs are colored based on the end position of the splice junctions, and Y-axis represents actual read counts.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    junctions = {}
    
    chrom, region_range = region.split(":")
    region_start, region_end = map(int, region_range.split("-"))

    # Process each BAM file
    for bam_file in bam_files:
        samfile = pysam.AlignmentFile(bam_file, "rb")
        for read in samfile.fetch(chrom, region_start, region_end):
            if read.is_unmapped or read.is_secondary:
                continue

            # Process splice junctions (CIGAR N operation)
            for cigar in read.cigartuples:
                if cigar[0] == 3:  # CIGAR N operation: skip region (splice junction)
                    start = read.reference_start
                    end = start + cigar[1]
                    if (start, end) not in junctions:
                        junctions[(start, end)] = 0
                    junctions[(start, end)] += 1

    # Filter junctions
    filtered_junctions = {(start, end): count for (start, end), count in junctions.items() if 18 <= start <= 100 and end <= 600}
    
    # Summary of junction counts
    total_junctions = len(junctions)
    print(f"Total Junctions: {total_junctions}")
    print(f"Filtered Junctions (Start between 18th and 100th bp): {len(filtered_junctions)}")

    # Write filtered output to a file
    output_file = "filtered_junctions_by_end.txt"
    with open("filtered_junction_output.txt", "w") as f:
        f.write("Start\tEnd\tRead Count\n")
        for (start, end), count in filtered_junctions.items():
            f.write(f"{start}\t{end}\t{count}\n")
    
    print(f"Filtered junctions saved to {output_file}")

    # Define a colormap based on the end positions of the junctions
    end_positions = [end for _, end in filtered_junctions.keys()]
    norm = mcolors.Normalize(vmin=min(end_positions), vmax=max(end_positions))
    cmap = cm.get_cmap("coolwarm")  # You can change the colormap if desired

    # Function to draw arcs with actual read counts
    def add_arc_with_color(start, end, count, color):
        arc_x = np.linspace(start, end, 100)
        # Use the actual read count for Y-axis height
        arc_y = count * np.sin(np.pi * (arc_x - start) / (end - start))  # Height based on actual read count
        ax.plot(arc_x, arc_y, color=color, lw=2)

    # Plot the arcs with colors based on end position
    for (start, end), count in filtered_junctions.items():
        color = cmap(norm(end))  # Get color based on the normalized end position
        add_arc_with_color(start, end, count, color)

    # Overlay exon structure from GTF
    exon_colors = {
        (1, 53): 'green',        # Exon 1 (Ckm)
        (54, 59): 'blue',        # Exon 2 (Intermediary)
        (60, 1266): 'purple'     # Exon 3 (F9)
    }
    
    exons = parse_gtf(gtf_file, region)
    for exon_start, exon_end in exons:
        color = exon_colors.get((exon_start, exon_end), 'gray')
        ax.add_patch(plt.Rectangle((exon_start, -0.2), exon_end - exon_start, 0.4, color=color, alpha=0.5))

    ax.set_xlabel(f"Genomic Position: {region}")
    ax.set_ylabel("Read Count (Junction Support)")  # Change label to reflect actual read counts
    ax.set_title(f"Sashimi Plot for Region {region} (Colored by End Position, Read Counts on Y-axis)")
    
    # Save plot as SVG
    plt.savefig("treatedRNAseq.svg", format="svg")
    plt.show()

# Usage example
bam_files = [
    "/hdd/InvivoRACE__concatenated_fastq/Invivo_RACE_alignment_redo/barcode08_sorted.bam",
    "/hdd/InvivoRACE__concatenated_fastq/Invivo_RACE_alignment_redo/barcode09_sorted.bam",
    "/hdd/InvivoRACE__concatenated_fastq/Invivo_RACE_alignment_redo/barcode10_sorted.bam",
    "/hdd/InvivoRACE__concatenated_fastq/Invivo_RACE_alignment_redo/barcode11_sorted.bam"
]

gtf_file = "/home/cnelsonlab/F9mRNARACE.gtf"
region = "F9mRNARACE:1-500"  # Define the region of interest
output_file = "/hdd/InvivoRACE__concatenated_fastq/Invivo_RACE_alignment_redo/RACE"

plot_sashimi_with_actual_read_counts(bam_files, gtf_file, region, output_file)


#####################shortread splice junction################
import pysam
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors

def parse_gtf(gtf_file, region):
    """
    Parse the GTF file to extract exon information for the region of interest.
    """
    exons = []
    chrom, region_range = region.split(":")
    region_start, region_end = map(int, region_range.split("-"))
    
    gtf = pd.read_csv(gtf_file, sep="\t", comment="#", header=None, names=[
        "seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
    
    for _, row in gtf.iterrows():
        if row['seqname'] == chrom and row['feature'] == "exon":
            if (row['start'] >= region_start and row['end'] <= region_end):
                exons.append((row['start'], row['end']))
    
    return exons

def plot_sashimi_based_on_exon_alignment(bam_files, gtf_file, region, output_file, min_exon1_length=20):
    """
    Generate a Sashimi plot based on alignment that covers at least 20 bp from Exon 1
    and starts an arc from Exon 1 end to Exon 2 or 3 alignment start.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    junctions = {}

    chrom, region_range = region.split(":")
    region_start, region_end = map(int, region_range.split("-"))

    exons = parse_gtf(gtf_file, region)

    # Extract Exon 1, Exon 2, and Exon 3 positions
    exon1_start, exon1_end = 1, 53  # Exon 1 (Ckm)
    exon2_start, exon2_end = 54, 59  # Exon 2 (Intermediary)
    exon3_start, exon3_end = 60, 1266  # Exon 3 (F9)

    # Process each BAM file
    for bam_file in bam_files:
        samfile = pysam.AlignmentFile(bam_file, "rb")
        for read in samfile.fetch(chrom, region_start, region_end):
            if read.is_unmapped or read.is_secondary:
                continue

            # Get read alignment details
            start = read.reference_start
            end = read.reference_end

            # Check if read covers at least 20 bp of Exon 1 (<=53)
            if start <= 53 and (min(53, end) - start >= min_exon1_length) and end > 53:
                
                # Now find where the alignment enters Exon 2 or Exon 3
                exon2_start_match = None
                exon3_start_match = None

                # Use CIGAR tuples to locate exon alignment precisely
                read_pos = start
                for cigar_op, cigar_len in read.cigartuples:
                    if cigar_op == 0:  # Match or mismatch
                        # Check if part of the read aligns to Exon 2
                        if exon2_start <= read_pos < exon2_end:
                            exon2_start_match = read_pos
                            break
                        # Check if part of the read aligns to Exon 3
                        elif exon3_start <= read_pos < exon3_end:
                            exon3_start_match = read_pos
                            break
                    read_pos += cigar_len  # Move read position based on CIGAR operation length

                # Determine the end of the arc (start of Exon 2 or 3)
                if exon2_start_match:
                    arc_end = exon2_start_match
                elif exon3_start_match:
                    arc_end = exon3_start_match
                else:
                    continue  # No valid match for Exon 2 or 3, skip

                # Draw the arc from end of Exon 1 alignment to start of Exon 2 or 3
                arc_start = min(53, end)  # Arc starts where alignment in Exon 1 ends
                if (arc_start, arc_end) not in junctions:
                    junctions[(arc_start, arc_end)] = 0
                junctions[(arc_start, arc_end)] += 1

    # Summary of junction counts
    total_junctions = len(junctions)
    print(f"Total Junctions: {total_junctions}")
    print(f"Filtered Junctions (with at least 20 bp in Exon 1 and matching Exon 2 or 3): {total_junctions}")

    # Handle the case where no junctions were found
    if not junctions:
        print("No valid junctions found with the specified criteria.")
        return

    # Write filtered output to a file
    output_file = "filtered_junctions_exon1_boundary.txt"
    with open(output_file, "w") as f:
        f.write("Start\tEnd\tRead Count\n")
        for (start, end), count in junctions.items():
            f.write(f"{start}\t{end}\t{count}\n")
    
    print(f"Filtered junctions saved to {output_file}")

    # Define a colormap based on the end positions of the junctions
    end_positions = [end for _, end in junctions.keys()]
    norm = mcolors.Normalize(vmin=min(end_positions), vmax=max(end_positions))
    cmap = cm.get_cmap("coolwarm")

    # Function to draw arcs with actual read counts
    def add_arc_with_color(start, end, count, color):
        arc_x = np.linspace(start, end, 100)
        arc_y = count * np.sin(np.pi * (arc_x - start) / (end - start))  # Height based on actual read count
        ax.plot(arc_x, arc_y, color=color, lw=2)

    # Plot the arcs with colors based on end position
    for (start, end), count in junctions.items():
        color = cmap(norm(end))  # Get color based on the normalized end position
        add_arc_with_color(start, end, count, color)

    # Overlay exon structure from GTF
    exon_colors = {
        (1, 53): 'green',        # Exon 1 (Ckm)
        (54, 59): 'blue',        # Exon 2 (Intermediary)
        (60, 1266): 'purple'     # Exon 3 (F9)
    }
    
    exons = parse_gtf(gtf_file, region)
    for exon_start, exon_end in exons:
        color = exon_colors.get((exon_start, exon_end), 'gray')
        ax.add_patch(plt.Rectangle((exon_start, -0.2), exon_end - exon_start, 0.4, color=color, alpha=0.5))

    ax.set_xlabel(f"Genomic Position: {region}")
    ax.set_ylabel("Read Count (Junction Support)")
    ax.set_title(f"Sashimi Plot Based on Exon Alignment for Region {region}")
    
    # Save plot as SVG
    plt.savefig("filtered_exon1_boundary_sashimi_adjusted.svg", format="svg")
    plt.show()


#bam_files = [
#    "/hdd/STAR/Output/Sample_5_Aligned.sortedByCoord.out.bam",
#    "/hdd/STAR/Output/Sample_6_Aligned.sortedByCoord.out.bam",
#    "/hdd/STAR/Output/Sample_7_Aligned.sortedByCoord.out.bam",
#    "/hdd/STAR/Output/Sample_8_Aligned.sortedByCoord.out.bam",
#    "/hdd/STAR/Output/Sample_13_Aligned.sortedByCoord.out.bam",
#    "/hdd/STAR/Output/Sample_14_Aligned.sortedByCoord.out.bam",
#    "/hdd/STAR/Output/Sample_15_Aligned.sortedByCoord.out.bam",
#    "/hdd/STAR/Output/Sample_16_Aligned.sortedByCoord.out.bam",
#]

#bam_files = [
#    "/hdd/InvivoRACE__concatenated_fastq/Invivo_RACE_alignment_redo/barcode08_sorted.bam",
#    "/hdd/InvivoRACE__concatenated_fastq/Invivo_RACE_alignment_redo/barcode09_sorted.bam",
#    "/hdd/InvivoRACE__concatenated_fastq/Invivo_RACE_alignment_redo/barcode10_sorted.bam",
#    "/hdd/InvivoRACE__concatenated_fastq/Invivo_RACE_alignment_redo/barcode11_sorted.bam"
#]


gtf_file = "/home/cnelsonlab/F9mRNARACE.gtf"
region = "F9mRNARACE:1-1266"  # Define the region of interest
output_file = "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/filtered_exon1_to_exon2_or_3_boundary_sashimi.pdf"

plot_sashimi_based_on_exon_alignment(bam_files, gtf_file, region, output_file)


##################manhattan visualizing structural variant########
import pysam
import vcfpy
import matplotlib.pyplot as plt

# Load the VCF file with structural variants
vcf_file = "/hdd/InvivoRACE__concatenated_fastq/Invivo_RACE_alignment_redo/min_5_combined_output_structural_variants.vcf"
vcf_reader = vcfpy.Reader.from_path(vcf_file)

# Store detected SVs (insertions, deletions, etc.)
sv_list = []
sv_types = {"DEL": [], "INS": [], "INV": [], "DUP": []}  # Separate storage for different types of SVs

for record in vcf_reader:
    if "SVTYPE" in record.INFO:
        sv_type = record.INFO['SVTYPE']
        chrom = record.CHROM
        pos = record.POS
        end = record.INFO['END'] if 'END' in record.INFO else pos
        
        # Add to corresponding SV type list
        if sv_type in sv_types:
            sv_types[sv_type].append((chrom, pos, end))

# Prepare Manhattan plot data (combine all SVs into one list)
manhattan_data = []
for sv_type, svs in sv_types.items():
    for sv in svs:
        chrom, start, end = sv
        sv_size = end - start
        manhattan_data.append((start, sv_size, sv_type))

# Plot Manhattan-style plot
plt.figure(figsize=(12, 6))

# Define colors for each type of SV
colors = {'DEL': 'red', 'INS': 'green', 'INV': 'blue', 'DUP': 'purple'}

# Plot each SV with its size
legend_added = set()
for start, sv_size, sv_type in manhattan_data:
    plt.scatter(start, sv_size, color=colors[sv_type], label=sv_type if sv_type not in legend_added else "", alpha=0.6, s=40)
    legend_added.add(sv_type)

# Annotate larger variants
for start, sv_size, sv_type in manhattan_data:
    if sv_size > 500:  # Only annotate large SVs
        plt.text(start, sv_size + 20, f'{sv_size} bp', fontsize=8, color='black')

# Set axis labels and title
plt.xlabel('Genomic Position')
plt.ylabel('SV Size (bp)')
plt.title('Manhattan Plot of Structural Variants by Size')

# Add legend (only one of each SV type)
plt.legend(loc='upper right', markerscale=2)

# Show plot
plt.tight_layout()
plt.savefig("vcf_RACE.svg", format = "svg")
plt.show()



#######combining sashimi output to manhattan code#####
import pysam
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import vcfpy
import time  # Import time module for profiling

def parse_gtf(gtf_file, region):
    exons = []
    chrom, region_range = region.split(":")
    region_start, region_end = map(int, region_range.split("-"))
    
    gtf = pd.read_csv(gtf_file, sep="\t", comment="#", header=None, names=[
        "seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
    
    for _, row in gtf.iterrows():
        if row['seqname'] == chrom and row['feature'] == "exon":
            if (row['start'] >= region_start and row['end'] <= region_end):
                exons.append((row['start'], row['end']))
    
    return exons

def parse_vcf(vcf_file):
    sv_list = []
    sv_reader = vcfpy.Reader.from_path(vcf_file)
    for record in sv_reader:
        if "SVTYPE" in record.INFO:
            sv_type = record.INFO['SVTYPE']
            chrom = record.CHROM
            pos = record.POS
            end = record.INFO['END'] if 'END' in record.INFO else pos
            sv_list.append((chrom, pos, end, sv_type))
    return sv_list

def plot_manhattan_sv_only(bam_files, gtf_file, vcf_file, region, output_file, min_exon1_length=20):
    fig, manhattan_ax = plt.subplots(figsize=(10, 6))  # Only Manhattan plot
    
    junctions = {}
    sv_overlaps = []

    chrom, region_range = region.split(":")
    region_start, region_end = map(int, region_range.split("-"))

    # Time parsing of GTF and VCF files
    start_time = time.time()
    exons = parse_gtf(gtf_file, region)
    print(f"Time for parsing GTF: {time.time() - start_time:.2f} seconds")

    start_time = time.time()
    sv_list = parse_vcf(vcf_file)  # Parse the VCF file for SVs
    print(f"Time for parsing VCF: {time.time() - start_time:.2f} seconds")

    # Process each BAM file
    for bam_file in bam_files:
        print(f"Processing BAM file: {bam_file}")
        samfile = pysam.AlignmentFile(bam_file, "rb")

        # Profile time for fetching reads
        start_time = time.time()
        read_count = 0
        for read in samfile.fetch(chrom, region_start, region_end):
            read_count += 1
            if read.is_unmapped or read.is_secondary:
                continue

            start = read.reference_start
            end = read.reference_end

            if start <= 53 and (min(53, end) - start >= min_exon1_length) and end > 53:
                # Perform CIGAR operation check (time-consuming)
                read_pos = start
                for cigar_op, cigar_len in read.cigartuples:
                    if cigar_op == 0:  # Match or mismatch
                        # Check alignment to Exon 2 or Exon 3
                        if 54 <= read_pos < 59 or 60 <= read_pos < 1266:
                            break
                    read_pos += cigar_len
                
                # Overlap check with SVs
                for sv in sv_list:
                    sv_chrom, sv_start, sv_end, sv_type = sv
                    if sv_chrom == chrom and (sv_start <= end and sv_end >= start):  # Overlap check
                        sv_overlaps.append((sv_start, sv_end, sv_type))

        print(f"Processed {read_count} reads from BAM file in {time.time() - start_time:.2f} seconds")

    # Plot Structural Variants (Manhattan plot)
    print(f"Plotting {len(sv_overlaps)} structural variants")
    start_time = time.time()
    
    colors = {'DEL': 'red', 'INS': 'green', 'INV': 'blue', 'DUP': 'purple'}
    for sv_start, sv_end, sv_type in sv_overlaps:
        sv_size = sv_end - sv_start
        manhattan_ax.scatter(sv_start, sv_size, color=colors.get(sv_type, 'black'), s=50)

    manhattan_ax.set_title('Manhattan Plot of Structural Variants Overlapping Sashimi Reads')
    manhattan_ax.set_xlabel('Genomic Position')
    manhattan_ax.set_ylabel('SV Size (bp)')

    plt.tight_layout()
    plt.savefig(output_file, format="pdf")
    plt.show()

    print(f"Plotting completed in {time.time() - start_time:.2f} seconds")

# Example usage with your BAM and VCF files
bam_files = [
    "/hdd/InvivoRACE__concatenated_fastq/Invivo_RACE_alignment_redo/combined_sorted.bam"]

gtf_file = "/home/cnelsonlab/F9mRNARACE.gtf"
vcf_file = "/hdd/InvivoRACE__concatenated_fastq/Invivo_RACE_alignment_redo/min_5_combined_output_structural_variants.vcf"
region = "F9mRNARACE:1-1266"  # Define the region of interest
output_file = "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/combined_sv_manhattan.pdf"

plot_manhattan_sv_only(bam_files, gtf_file, vcf_file, region, output_file)


####upside down sashimi plot####
import pysam
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors

def parse_gtf(gtf_file, region):
    """
    Parse the GTF file to extract exon information for the region of interest.
    """
    exons = []
    chrom, region_range = region.split(":")
    region_start, region_end = map(int, region_range.split("-"))
    
    gtf = pd.read_csv(gtf_file, sep="\t", comment="#", header=None, names=[
        "seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
    
    for _, row in gtf.iterrows():
        if row['seqname'] == chrom and row['feature'] == "exon":
            if (row['start'] >= region_start and row['end'] <= region_end):
                exons.append((row['start'], row['end']))
    
    return exons

def plot_sashimi_based_on_exon_alignment(bam_files_long, bam_files_short, gtf_file, region, output_file, min_exon1_length=20):
    """
    Generate a Sashimi plot with long-read data on top and short-read data upside down on the bottom
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), sharex=True)

    chrom, region_range = region.split(":")
    region_start, region_end = map(int, region_range.split("-"))

    exons = parse_gtf(gtf_file, region)

    # Extract Exon 1, Exon 2, and Exon 3 positions
    exon1_start, exon1_end = 1, 53  # Exon 1 (Ckm)
    exon2_start, exon2_end = 54, 59  # Exon 2 (Intermediary)
    exon3_start, exon3_end = 60, 1266  # Exon 3 (F9)

    def process_bam_files(bam_files, ax, inverted=False):
        junctions = {}

        # Process each BAM file
        for bam_file in bam_files:
            samfile = pysam.AlignmentFile(bam_file, "rb")
            for read in samfile.fetch(chrom, region_start, region_end):
                if read.is_unmapped or read.is_secondary:
                    continue

                # Get read alignment details
                start = read.reference_start
                end = read.reference_end

                # Check if read covers at least 20 bp of Exon 1 (<=53)
                if start <= 53 and (min(53, end) - start >= min_exon1_length) and end > 53:
                    
                    # Now find where the alignment enters Exon 2 or Exon 3
                    exon2_start_match = None
                    exon3_start_match = None

                    # Use CIGAR tuples to locate exon alignment precisely
                    read_pos = start
                    for cigar_op, cigar_len in read.cigartuples:
                        if cigar_op == 0:  # Match or mismatch
                            # Check if part of the read aligns to Exon 2
                            if exon2_start <= read_pos < exon2_end:
                                exon2_start_match = read_pos
                                break
                            # Check if part of the read aligns to Exon 3
                            elif exon3_start <= read_pos < exon3_end:
                                exon3_start_match = read_pos
                                break
                        read_pos += cigar_len  # Move read position based on CIGAR operation length

                    # Determine the end of the arc (start of Exon 2 or 3)
                    if exon2_start_match:
                        arc_end = exon2_start_match
                    elif exon3_start_match:
                        arc_end = exon3_start_match
                    else:
                        continue  # No valid match for Exon 2 or 3, skip

                    # Draw the arc from end of Exon 1 alignment to start of Exon 2 or 3
                    arc_start = min(53, end)  # Arc starts where alignment in Exon 1 ends
                    if (arc_start, arc_end) not in junctions:
                        junctions[(arc_start, arc_end)] = 0
                    junctions[(arc_start, arc_end)] += 1

        # Plot arcs with colors based on end position
        end_positions = [end for _, end in junctions.keys()]
        if end_positions:
            norm = mcolors.Normalize(vmin=min(end_positions), vmax=max(end_positions))
            cmap = cm.get_cmap("coolwarm")

            # Function to draw arcs with actual read counts
            def add_arc_with_color(start, end, count, color):
                arc_x = np.linspace(start, end, 100)
                arc_y = count * np.sin(np.pi * (arc_x - start) / (end - start))  # Height based on actual read count
                if inverted:
                    arc_y = -arc_y  # Invert arc for the second plot
                ax.plot(arc_x, arc_y, color=color, lw=2)

            # Plot the arcs with colors based on end position
            for (start, end), count in junctions.items():
                color = cmap(norm(end))  # Get color based on the normalized end position
                add_arc_with_color(start, end, count, color)

        return junctions

    # Process long-read BAM files (top plot)
    process_bam_files(bam_files_long, ax1)

    # Process short-read BAM files (bottom plot, inverted)
    process_bam_files(bam_files_short, ax2, inverted=True)

    # Overlay exon structure from GTF
    exon_colors = {
        (1, 53): 'green',        # Exon 1 (Ckm)
        (54, 59): 'blue',        # Exon 2 (Intermediary)
        (60, 1266): 'purple'     # Exon 3 (F9)
    }

    for exon_start, exon_end in exons:
        color = exon_colors.get((exon_start, exon_end), 'gray')
        ax1.add_patch(plt.Rectangle((exon_start, -0.2), exon_end - exon_start, 0.4, color=color, alpha=0.5))
        ax2.add_patch(plt.Rectangle((exon_start, -0.2), exon_end - exon_start, 0.4, color=color, alpha=0.5))

    ax1.set_title("Long-read Sashimi Plot")
    ax2.set_title("Short-read Sashimi Plot (Upside Down)")

    ax2.set_xlabel(f"Genomic Position: {region}")
    ax1.set_ylabel("Read Count (Junction Support)")
    ax2.set_ylabel("Read Count (Inverted)")

    plt.savefig(output_file, format="svg")
    plt.show()

# Usage example
bam_files_long = [
    "/hdd/InvivoRACE__concatenated_fastq/Invivo_RACE_alignment_redo/barcode08_sorted.bam",
    "/hdd/InvivoRACE__concatenated_fastq/Invivo_RACE_alignment_redo/barcode09_sorted.bam",
]

bam_files_short = [
    "/hdd/STAR/Output/Sample_2_Aligned.sortedByCoord.out.bam",
    "/hdd/STAR/Output/Sample_3_Aligned.sortedByCoord.out.bam",
]

gtf_file = "/home/cnelsonlab/F9mRNARACE.gtf"
region = "F9mRNARACE:1-1266"  # Define the region of interest
output_file = "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/long_vs_short_sashimi.svg"

plot_sashimi_based_on_exon_alignment(bam_files_long, bam_files_short, gtf_file, region, output_file)
