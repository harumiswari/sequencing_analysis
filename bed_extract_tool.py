import pandas as pd
import subprocess
import os

def convert_csv_to_bed(csv_filepath, bed_filepath, chrom_map):
    # Load CSV data
    csv_data = pd.read_csv(csv_filepath)

    # Map chromosome names using provided dictionary
    csv_data['Chrom'] = csv_data['Chrom'].map(chrom_map)

    # Calculate start and end positions for BED format
    csv_data['Start'] = csv_data['Position'] - 31  # Extend 30 bp upstream, adjust to 0-based start
    csv_data['End'] = csv_data['Position'] + 30    # Extend 30 bp downstream

    # Select relevant columns for the BED file
    bed_data = csv_data[['Chrom', 'Start', 'End', 'Counts']]
    
    # Save to BED file
    bed_data.to_csv(bed_filepath, sep='\t', index=False, header=False)
    
    # Display the first few rows to confirm
    print(bed_data.head())

def run_bedtools_intersect(bed_file, gtf_file, output_file):
    # Check if the BED file exists and is not empty
    if os.path.exists(bed_file) and os.path.getsize(bed_file) > 0:
        bedtools_command = [
            'bedtools', 'intersect',
            '-a', bed_file,
            '-b', gtf_file,
            '-wa', '-wb'
        ]
        
        try:
            result = subprocess.run(bedtools_command, capture_output=True, text=True, check=True)
            with open(output_file, 'w') as file:
                file.write(result.stdout)
            print(f"Output written to {output_file}")
        except subprocess.CalledProcessError as e:
            print("Error in bedtools intersect:", e.stderr)
    else:
        print(f"Failed to proceed: {bed_file} does not exist or is empty.")

# Chromosome mapping dictionary
chrom_map = {
    'CM000994.3': '1', 'CM000995.3': '2', 'CM000996.3': '3', 'CM000997.3': '4',
    'CM000998.3': '5', 'CM000999.3': '6', 'CM001000.3': '7', 'CM001001.3': '8',
    'CM001002.3': '9', 'CM001003.3': '10', 'CM001004.3': '11', 'CM001005.3': '12',
    'CM001006.3': '13', 'CM001007.3': '14', 'CM001008.3': '15', 'CM001009.3': '16',
    'CM001010.3': '17', 'CM001011.3': '18', 'CM001012.3': '19', 'CM001013.3': 'X',
    'CM001014.3': 'Y'
}

gtf_file = '/Users/gwisna/Mus_musculus.GRCm39.112.gtf'

# Process multiple barcodes
for barcode in range(31, 41): 
    csv_filepath = f'/Users/gwisna/Desktop/hF9mousecombined/BC{barcode}/BC{barcode}_integration_site_distribution_filtered.csv'
    bed_filepath = f'/Users/gwisna/Desktop/hF9mousecombined/BC{barcode}/BC{barcode}_updated_positions.bed'
    output_file = f'/Users/gwisna/Desktop/hF9mousecombined/BC{barcode}/BC{barcode}_annotated.txt'
    
    # Convert CSV to BED
    convert_csv_to_bed(csv_filepath, bed_filepath, chrom_map)
    
    # Run bedtools intersect
    run_bedtools_intersect(bed_filepath, gtf_file, output_file)


######To calculate the hits########

def calculate_hits(file_path, output_file):
    # Load data
    data = pd.read_csv(file_path, sep='\t', header=None)
    
    # Count total hits per chromosome
    total_hits = data.groupby(0)[3].sum().reset_index()
    total_hits.columns = ['Chromosome', 'Total_Hits']
    
    # Filter for 'gene' in the 7th column and count hits
    gene_hits = data[data[6].str.contains('gene', case=False, na=False)]
    gene_specific_hits = gene_hits.groupby(0)[3].sum().reset_index()
    gene_specific_hits.columns = ['Chromosome', 'Gene_Hits']

    # Merge the total hits and gene-specific hits into one dataframe
    merged_data = pd.merge(total_hits, gene_specific_hits, on='Chromosome', how='left')
    merged_data['Gene_Hits'] = merged_data['Gene_Hits'].fillna(0)  # Replace NaN with 0 where there are no gene hits

    # Calculate Final Hits by subtracting Gene Hits from Total Hits
    merged_data['Final_Hits'] = merged_data['Total_Hits'] - merged_data['Gene_Hits']

    # Save the results to file
    merged_data.to_csv(output_file, index=False)
    print(f"Output written to {output_file}")

# Process multiple barcodes
for barcode in range(31, 41):  # Adjust the range as needed
    input_file = f'/Users/gwisna/Desktop/hF9mousecombined/BC{barcode}/BC{barcode}_annotated.txt'
    output_file = f'/Users/gwisna/Desktop/hF9mousecombined/BC{barcode}/BC{barcode}_calculated1.csv'
    calculate_hits(input_file, output_file)
