#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 15:44:05 2024

@author: cnelsonlab
"""

import pysam
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

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

def plot_sashimi(bam_files, gtf_file, region, output_file):
    """
    Generate a Sashimi plot from multiple BAM files and a GTF annotation file.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    junctions = {}
    insertions = {}  # Track insertions
    deletions = {}   # Track deletions

    chrom, region_range = region.split(":")
    region_start, region_end = map(int, region_range.split("-"))

    # Process each BAM file
    for bam_file in bam_files:
        samfile = pysam.AlignmentFile(bam_file, "rb")
        for read in samfile.fetch(chrom, region_start, region_end):
            if read.is_unmapped or read.is_secondary:
                continue

            # Skip reads shorter than 500 bp
            #if read.query_length < 500:
            #    continue

            # Skip reads with soft clipping (CIGAR 'S')
            #if any(cigar[0] == 4 for cigar in read.cigartuples):  # CIGAR S operation: soft clipping
            #    continue

            # Process splice junctions (CIGAR N operation) and track indels
            for cigar in read.cigartuples:
                if cigar[0] == 3:  # CIGAR N operation: skip region (splice junction)
                    start = read.reference_start
                    end = start + cigar[1]
                    if (start, end) not in junctions:
                        junctions[(start, end)] = 0
                    junctions[(start, end)] += 1
                elif cigar[0] == 1:  # CIGAR I operation: insertion
                    insertion_pos = read.reference_start
                    if insertion_pos not in insertions:
                        insertions[insertion_pos] = 0
                    insertions[insertion_pos] += 1
                elif cigar[0] == 2:  # CIGAR D operation: deletion
                    deletion_pos = read.reference_start
                    if deletion_pos not in deletions:
                        deletions[deletion_pos] = 0
                    deletions[deletion_pos] += 1

    # Filter junctions with more than 5 supporting reads
    #filtered_junctions = {key: count for key, count in junctions.items() if start == 10}
    filtered_junctions = {(start, end): count for (start, end), count in junctions.items() if 18 <= start <= 100 and end <= 600}
   
    # Summary of junction counts
    total_junctions = len(junctions)
    filtered_junctions_count = len(filtered_junctions)
    
    print(f"Total Junctions: {total_junctions}")
    #print(f"Filtered Junctions (Count > 10): {filtered_junctions_count}")
    print(f"Filtered Junctions (Start between 18th and 100th bp): {filtered_junctions_count}")

    # Write full output to a file
    output_file = "filtered_junctions.txt"
    with open("filtered_junction_output.txt", "w") as f:
    #    f.write("All Junctions with Counts:\n")
    #    for (start, end), count in junctions.items():
    #        f.write(f"Junction {start}-{end}, Count: {count}\n")
        
        f.write("Start\tEnd\tRead Count\n")
        for (start, end), count in filtered_junctions.items():
            f.write(f"{start}\t{end}\t{count}\n")
    print(f"Filtered junctions saved to {output_file}")

    # Function to draw parabolic arcs for the splice junctions
    def add_arc(start, end, count):
        arc_x = np.linspace(start, end, 100)
        arc_y = count - (count * (arc_x - start) * (arc_x - end) / ((end - start)**2))  # Parabolic arc
        ax.plot(arc_x, arc_y, color="blue", lw=2)

    # Plot the arcs for each splice junction from filtered junctions
    for (start, end), count in filtered_junctions.items():  # Now using 'filtered_junctions'
        add_arc(start, end, count)

    # Plot insertions as red vertical lines
    #for pos, count in insertions.items():
    #    ax.axvline(pos, color="red", linestyle="--", lw=1, alpha=0.5, label="Insertion")

    # Plot deletions as orange vertical lines
    #for pos, count in deletions.items():
    #    ax.axvline(pos, color="orange", linestyle="--", lw=1, alpha=0.5, label="Deletion")

    exon_colors = {
    (1, 53): 'green',        # Exon 1 (Ckm)
    (54, 59): 'blue',        # Exon 2 (Intermediary)
    (60, 1266): 'purple'     # Exon 3 (F9)
}

    # Overlay exon structure from GTF
    exons = parse_gtf(gtf_file, region)
    for exon_start, exon_end in exons:
        color = exon_colors.get((exon_start, exon_end), 'gray')
        ax.add_patch(plt.Rectangle((exon_start, 0), exon_end - exon_start, 1, color=color, alpha=0.5))

    ax.set_xlabel(f"Genomic Position: {region}")
    ax.set_ylabel("Read Count (Junction Support)")
    ax.set_title(f"Sashimi Plot for Region {region}")
    
    # Save plot to file
    plt.savefig("treatedRNAseq.svg", format="svg")
    plt.show()

# Usage example
bam_files = [
    "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/Sample_1.bam",
    "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/Sample_2.bam",
    "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/Sample_3.bam",
    "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/Sample_4.bam"
]

gtf_file = "/home/cnelsonlab/F9mRNARACE.gtf"
region = "F9mRNARACE:1-1266"  # Define the region of interest
output_file = "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/treatedRNAseq.pdf"

plot_sashimi(bam_files, gtf_file, region, output_file)

gtf = pd.read_csv(gtf_file, sep="\t", comment="#", header=None, names=[
    "seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
print(gtf.head())  # This will print the first few rows of your GTF file

exons = parse_gtf(gtf_file, region)
print("Parsed exons:", exons)


#############for multi track#############

import pysam
import matplotlib.pyplot as plt
import numpy as np

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

def plot_combined_coverage_and_junctions(bam_files, gtf_file, region, output_file):
    """
    Generate a combined plot of coverage and splice junctions across multiple samples.
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

bam_files = [
    "/hdd/InvivoRACE__concatenated_fastq/Invivo_RACE_alignment_redo/barcode08_sorted.bam",
    "/hdd/InvivoRACE__concatenated_fastq/Invivo_RACE_alignment_redo/barcode09_sorted.bam",
    "/hdd/InvivoRACE__concatenated_fastq/Invivo_RACE_alignment_redo/barcode10_sorted.bam",
    "/hdd/InvivoRACE__concatenated_fastq/Invivo_RACE_alignment_redo/barcode11_sorted.bam"
]

gtf_file = "/home/cnelsonlab/F9mRNARACE.gtf"
region = "F9mRNARACE:1-500"  # Define the region of interest
output_file = "multiple_tract_coverage.svg"

plot_combined_coverage_and_junctions(bam_files, gtf_file, region, output_file)


################combining arc into 1 plot###################
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
    "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/Sample_1.bam",
    "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/Sample_2.bam",
    "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/Sample_3.bam",
    "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/Sample_4.bam"
]

gtf_file = "/home/cnelsonlab/F9mRNARACE.gtf"
region = "F9mRNARACE:1-500"  # Define the region of interest
output_file = "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/treatedRNAseq.pdf"

plot_sashimi_with_actual_read_counts(bam_files, gtf_file, region, output_file)


#####################shortread################
import pysam
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors

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

def plot_sashimi_based_on_exon_alignment(bam_files, gtf_file, region, output_file, min_exon1_length=20):
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

            # Get read alignment details
            start = read.reference_start
            end = read.reference_end

            # Check if read covers at least 20 bp of Exon 1 and extends past 53 bp
            if start <= 53 and (end - start >= min_exon1_length) and end > 53:
                # Consider this read as valid and create an arc
                if (start, end) not in junctions:
                    junctions[(start, end)] = 0
                junctions[(start, end)] += 1

    # Summary of junction counts
    total_junctions = len(junctions)
    print(f"Total Junctions: {total_junctions}")
    print(f"Filtered Junctions (with at least 20 bp in Exon 1 and crossing 53 bp): {total_junctions}")

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
    cmap = cm.get_cmap("coolwarm")  # You can change the colormap if desired

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
    ax.set_xticks(np.arange(0, region_end + 50, 50))
    # Save plot as SVG
    plt.savefig("filtered_exon1_boundary_sashimi.svg", format="svg")
    plt.show()

# Usage example
#bam_files = [
#    "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/Sample_1.bam",
#    "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/Sample_2.bam",
#    "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/Sample_3.bam",
#    "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/Sample_4.bam"
#]

bam_files = [
    "/hdd/STAR/Output/Sample_2_Aligned.sortedByCoord.out.bam",
    "/hdd/STAR/Output/Sample_3_Aligned.sortedByCoord.out.bam",
    "/hdd/STAR/Output/Sample_4_Aligned.sortedByCoord.out.bam"
]

gtf_file = "/home/cnelsonlab/F9mRNARACE.gtf"
region = "F9mRNARACE:1-1266"  # Define the region of interest
output_file = "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/filtered_exon1_boundary_sashimi.pdf"

plot_sashimi_based_on_exon_alignment(bam_files, gtf_file, region, output_file)



#######three output#############
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

def plot_sashimi_based_on_exon_alignment(bam_files, gtf_file, region, output_file, min_exon1_length=18):
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

# Usage example
#bam_files = [
#    "/hdd/STAR/Output/Sample_2_Aligned.sortedByCoord.out.bam",
#    "/hdd/STAR/Output/Sample_3_Aligned.sortedByCoord.out.bam",
#    "/hdd/STAR/Output/Sample_4_Aligned.sortedByCoord.out.bam"
#]

bam_files = [
    "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/Sample_1.bam",
    "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/Sample_2.bam",
    "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/Sample_3.bam",
    "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/Sample_4.bam"
]


gtf_file = "/home/cnelsonlab/F9mRNARACE.gtf"
region = "F9mRNARACE:1-1266"  # Define the region of interest
output_file = "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/filtered_exon1_to_exon2_or_3_boundary_sashimi.pdf"

plot_sashimi_based_on_exon_alignment(bam_files, gtf_file, region, output_file)






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

def plot_sashimi_based_on_exon_alignment(bam_files, gtf_file, region, output_file, readcount_output_file, min_exon1_length=20):
    """
    Generate a Sashimi plot based on alignment that covers at least 20 bp from Exon 1
    and starts an arc from Exon 1 end to Exon 2 or 3 alignment start. Also, write read count details to an output file.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    junctions = {}
    readcount_info = []  # List to store read count information

    chrom, region_range = region.split(":")
    region_start, region_end = map(int, region_range.split("-"))

    exons = parse_gtf(gtf_file, region)

    # Extract Exon 1, Exon 2, and Exon 3 positions (adjusted for 0-based indexing)
    exon1_start, exon1_end = 0, 52  # Exon 1 (Ckm) is now 0-based (1-53 -> 0-52)
    exon2_start, exon2_end = 53, 58  # Exon 2 (Intermediary) is now 0-based (54-59 -> 53-58)
    exon3_start, exon3_end = 59, 1265  # Exon 3 (F9) is now 0-based (60-1266 -> 59-1265)

    # Process each BAM file
    for bam_file in bam_files:
        samfile = pysam.AlignmentFile(bam_file, "rb")
        for read in samfile.fetch(chrom, region_start, region_end):
            if read.is_unmapped or read.is_secondary:
                continue

            # Get read alignment details
            aligned_positions = read.get_reference_positions()

            # Check if read covers at least 20 bp of Exon 1 (<=53)
            exon1_positions = [pos for pos in aligned_positions if exon1_start <= pos <= exon1_end]
            if len(exon1_positions) < min_exon1_length:
                continue
               
                # Now find where the alignment enters Exon 2 or Exon 3
                exon2_start_match = next((pos for pos in aligned_positions if exon2_start <= pos <= exon2_end), None)
                exon3_start_match = next((pos for pos in aligned_positions if exon3_start <= pos <= exon3_end), None)

 
                # Determine the end of the arc (start of Exon 2 or 3)
                if exon2_start_match or exon3_start_match:
                    arc_start = max(exon1_positions)
                    arc_end = exon2_start_match if exon2_start_match else exon3_start_match
                    if (arc_start, arc_end) not in junctions: 
                        junctions[(arc_start, arc_end)] = 0
                    junctions[(arc_start, arc_end)] += 1

                # Save read count information for Exon 1 and Exon 2/3 alignment
                readcount_info.append(f"Read ID: {read.query_name}, Exon 1 End: {arc_start}, Exon 2/3 Start: {arc_end}")

    # Summary of junction counts
    total_junctions = len(junctions)
    print(f"Total Junctions: {total_junctions}")
    print(f"Filtered Junctions (with at least 20 bp in Exon 1 and matching Exon 2 or 3): {total_junctions}")

    # Handle the case where no junctions were found
    if not junctions:
        print("No valid junctions found with the specified criteria.")
        return

    # Write filtered output (start and end of arcs) to a file
    with open(output_file, "w") as f:
        f.write("Start\tEnd\tRead Count\n")
        for (start, end), count in junctions.items():
            f.write(f"{start}\t{end}\t{count}\n")
    
    print(f"Filtered junctions saved to {output_file}")

    # Write read count information to another file
    with open(readcount_output_file, "w") as f:
        f.write("Reads aligning to Exon 1 and Exon 2/3:\n")
        for read_info in readcount_info:
            f.write(read_info + "\n")

    print(f"Read count information saved to {readcount_output_file}")

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

# Usage example
#bam_files = [
#    "/hdd/STAR/Output/Sample_2_Aligned.sortedByCoord.out.bam",
#    "/hdd/STAR/Output/Sample_3_Aligned.sortedByCoord.out.bam",
#    "/hdd/STAR/Output/Sample_4_Aligned.sortedByCoord.out.bam"
#]

bam_files = [
    "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/Sample_1.bam",
    "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/Sample_2.bam",
    "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/Sample_3.bam",
    "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/Sample_4.bam"
]

gtf_file = "/home/cnelsonlab/F9mRNARACE.gtf"
region = "F9mRNARACE:1-1266"  # Define the region of interest
output_file = "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/filtered_exon1_to_exon2_or_3_boundary_sashimi.txt"
readcount_output_file = "/hdd/DRS/Salmon_shortread/F9mRNARACE_aligned/exon1_to_exon2_or_3_readcount.txt"

plot_sashimi_based_on_exon_alignment(bam_files, gtf_file, region, output_file, readcount_output_file)


