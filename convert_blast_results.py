###python convert_blast_results.py QTL_A02.fasta.ECD10.7.out QTL_A02.fasta.ECD10.csv -chr A02 -e 0.0001 -b 20000
import argparse
import csv
import os

def convert_blast_results(input_file_path, output_file_path, chromosome_filter, evalue_threshold, bit_score_threshold, sort_by_start):
    # Check if the input file exists
    if not os.path.isfile(input_file_path):
        print(f"Error: The file {input_file_path} does not exist.")
        return
    
    # Read the input file
    with open(input_file_path, 'r') as file:
        lines = file.readlines()

    # Remove the first 4 lines (header lines) and the row containing field information
    filtered_lines = lines[4:]  # This removes the first 4 lines
    filtered_lines = [line for line in filtered_lines if "# Fields:" not in line]

    # Get the header from the row that starts with "# Fields:" and split by commas
    header_line = lines[3].strip()  # This is the line with the column titles
    header = header_line.strip("# ").split(", ")

    # Add the row with the hits count (line # 6)
    hits_line = lines[5].strip()
    filtered_lines.insert(0, hits_line)

    # Extract data lines (excluding the first line of hits information and the header row)
    data_lines = [line.strip() for line in filtered_lines if line and not line.startswith('#')]

    # Filter data lines based on the chromosome, evalue, and bit score criteria
    filtered_data = []
    for line in data_lines:
        columns = line.split()
        chromosome = columns[1]
        evalue = float(columns[10]) if columns[10] != "0.0" else 0.0  # Avoid issues with 0.0
        bit_score = float(columns[11])

        # Apply chromosome filter
        if chromosome_filter and chromosome != chromosome_filter:
            continue
        
        # Apply evalue filter
        if evalue_threshold and evalue > float(evalue_threshold):
            continue
        
        # Apply bit score filter
        if bit_score_threshold and bit_score < float(bit_score_threshold):
            continue
        
        filtered_data.append(columns)

    # Sort the filtered data by "s. start" (9th column) if required
    if sort_by_start:
        filtered_data.sort(key=lambda x: int(x[8]))  # s. start is in the 9th column (index 8)

    # Write the filtered data to a CSV file
    with open(output_file_path, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        # Write the header
        csvwriter.writerow(['# hits found'] + header)  # Include hits line as a column
        # Write the filtered data rows
        for row in filtered_data:
            csvwriter.writerow(row)

    print(f"Conversion completed. Output saved to {output_file_path}")

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Convert BLAST results to CSV format with filtering and sorting options.")
    parser.add_argument('infile', help="Path to the input BLAST results file.")
    parser.add_argument('outfile', help="Path to the output CSV file.")
    
    # Add optional arguments for filtering
    parser.add_argument('-chr', '--chromosome', type=str, help="Chromosome to filter by (e.g., A02).")
    parser.add_argument('-e', '--evalue', type=float, help="Maximum evalue threshold to filter by.")
    parser.add_argument('-b', '--bit_score', type=float, help="Minimum bit score threshold to filter by.")
    
    # Add argument for sorting
    parser.add_argument('-s', '--sort', type=str, choices=['true', 'false'], default='true', 
                        help="Whether to sort by s. start (default: true). Set to 'false' to disable sorting.")
    
    args = parser.parse_args()

    # Convert BLAST results with filters and sorting
    convert_blast_results(args.infile, args.outfile, args.chromosome, args.evalue, args.bit_score, args.sort == 'true')

if __name__ == "__main__":
    main()
