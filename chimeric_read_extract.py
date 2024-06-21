import pysam
import pandas as pd
import matplotlib.pyplot as plt

def categorize_chimeric_reads(samfile, chrom, start, end):
    filtered_reads = []
    ITRtoITRinsert_reads = []
    read_names_processed = set()  # Set to track processed read names

    for read in samfile.fetch(chrom, start, end):
        if read.query_name in read_names_processed:
            continue  # Skip this read if it has already been processed
        if read.has_tag('SA'):
            read_names_processed.add(read.query_name)
            original_chrom = read.reference_name
            original_pos = read.reference_start
            chimeric_info = read.get_tag('SA')
            for chimeric_part in chimeric_info.split(';'):
                if chimeric_part:
                    c_chrom, c_pos, c_strand, c_cigar, c_mapq, c_nm = chimeric_part.split(',')
                    if c_chrom != original_chrom:
                        filtered_reads.append((read.query_name, f"{original_chrom}:{original_pos}", f"{c_chrom}:{c_pos} {c_strand}", c_cigar))
                    if c_chrom == chrom:
                        ITRtoITRinsert_reads.append((read.query_name, f"{original_chrom}:{original_pos}", f"{c_chrom}:{c_pos} {c_strand}", c_cigar))
    return filtered_reads, ITRtoITRinsert_reads

def process_barcode(barcode):
    base_path = f'/Users/gwisna/Desktop/hF9mousecombined/BC{barcode}/'
    bam_path = f'{base_path}BC{barcode}_deduplicated.bam'
    samfile = pysam.AlignmentFile(bam_path, "rb")

    chromosome_name = 'ITRtoITRinsert'
    region_start = 600
    region_end = 900

    filtered_chimeric_reads, ITRtoITRinsert_reads = categorize_chimeric_reads(samfile, chromosome_name, region_start, region_end)

    # Save filtered chimeric reads
    with open(f'{base_path}BC{barcode}_chimeric_reads_results_filtered.txt', 'w') as file:
        file.write("Read Name\tOriginal Align\tChimeric Align\tCigar\n")
        for read in filtered_chimeric_reads:
            file.write(f"{read[0]}\t{read[1]}\t{read[2]}\t{read[3]}\n")

    # Save reads specific to 'ITRtoITRinsert'
    with open(f'{base_path}BC{barcode}_chimeric_reads_results_ITRtoITRinsert.txt', 'w') as file:
        file.write("Read Name\tOriginal Align\tChimeric Align\tCigar\n")
        for read in ITRtoITRinsert_reads:
            file.write(f"{read[0]}\t{read[1]}\t{read[2]}\t{read[3]}\n")

    samfile.close()

    # Process and analyze filtered data
    data_filtered = pd.read_csv(f'{base_path}BC{barcode}_chimeric_reads_results_filtered.txt', sep='\t')
    data_filtered['Priority'] = data_filtered['Chimeric Align'].apply(lambda x: 1 if 'CM001000.3' in x else 0)
    data_filtered.sort_values(by='Priority', ascending=False, inplace=True)
    data_filtered.drop_duplicates(subset='Read Name', keep='first', inplace=True)
    data_filtered.drop(columns=['Priority'], inplace=True)  # Clean up temporary column
    data_filtered['Chrom'] = data_filtered['Chimeric Align'].apply(lambda x: x.split(':')[0])
    data_filtered['Position'] = data_filtered['Chimeric Align'].apply(lambda x: int(x.split(':')[1].split(' ')[0]))
    location_counts_filtered = data_filtered.groupby(['Chrom', 'Position']).size().reset_index(name='Counts')
    location_counts_filtered.to_csv(f'{base_path}BC{barcode}_integration_site_distribution_filtered.csv', index=False)

    # Process and analyze ITRtoITRinsert data
    data_ITRtoITRinsert = pd.read_csv(f'{base_path}BC{barcode}_chimeric_reads_results_ITRtoITRinsert.txt', sep='\t')
    data_ITRtoITRinsert['Priority'] = data_ITRtoITRinsert['Chimeric Align'].apply(lambda x: 1 if 'CM001000.3' in x else 0)
    data_ITRtoITRinsert.sort_values(by='Priority', ascending=False, inplace=True)
    data_ITRtoITRinsert.drop_duplicates(subset='Read Name', keep='first', inplace=True)
    data_ITRtoITRinsert.drop(columns=['Priority'], inplace=True)  # Clean up temporary column
    data_ITRtoITRinsert['Chrom'] = data_ITRtoITRinsert['Chimeric Align'].apply(lambda x: x.split(':')[0])
    data_ITRtoITRinsert['Position'] = data_ITRtoITRinsert['Chimeric Align'].apply(lambda x: int(x.split(':')[1].split(' ')[0]))
    location_counts_ITRtoITRinsert = data_ITRtoITRinsert.groupby(['Chrom', 'Position']).size().reset_index(name='Counts')
    location_counts_ITRtoITRinsert.to_csv(f'{base_path}BC{barcode}_integration_site_distribution_ITRtoITRinsert.csv', index=False)

# Loop over bar


# Loop over barcodes from 31 to 40
for barcode in range(31, 41):
    process_barcode(barcode)

######To get fastq file from filtered bam#######

def filter_fastq(fastq_file, read_names, output_file):
    """Filter a FASTQ file based on a set of read names and write the matching records to a new file."""
    with open(fastq_file, 'r') as fq, open(output_file, 'w') as out_fq:
        while True:
            header = fq.readline()  # Read the header line
            if not header:
                break  # End of file reached
            sequence = fq.readline()  # Sequence line
            plus = fq.readline()  # '+' line
            quality = fq.readline()  # Quality line

            # Extract the read name from the header and check if it's in the provided set
            read_name_in_header = header.split()[0][1:]  # Assumes read name is the first part of the header, removes '@'
            if read_name_in_header in read_names:
                out_fq.write(header + sequence + plus + quality)

# Loop over barcodes from 31 to 40
for barcode in range(31, 41):
    # Load the read names from a CSV file where read names are stored in a column named 'Read Name'
    data_path = f'/Users/gwisna/Desktop/hF9mousecombined/BC{barcode}/BC{barcode}_chimeric_reads_results_ITRtoITRinsert.txt'
    data = pd.read_csv(data_path, delimiter='\t')
    read_names = set(data['Read Name'])  # Create a set of read names for efficient searching

    # Specify the path to the original FASTQ file and the output file
    fastq_file = f'/Users/gwisna/Desktop/with__pear/BC{barcode}/BC{barcode}_umi.fastq'
    output_file = f'/Users/gwisna/Desktop/hF9mousecombined/BC{barcode}/BC{barcode}_chimeric_reads_results_ITRtoITRinsert.fastq'

    # Filter the FASTQ file based on the read names extracted from the CSV
    filter_fastq(fastq_file, read_names, output_file)
