import argparse
from Bio import SeqIO

def parse_ids(id_file):
    """Reads the list file and returns a set of IDs to search for."""
    with open(id_file, 'r') as file:
        return set(line.strip() for line in file)

def extract_fasta(input_fasta, id_list, output_fasta):
    """Extracts sequences from input_fasta that match IDs in id_list and writes to output_fasta."""
    # Read the list of IDs to match
    ids_to_extract = parse_ids(id_list)
    
    # Open output file and write matching records
    with open(output_fasta, 'w') as output_handle:
        # Parse through each record in the FASTA file
        for record in SeqIO.parse(input_fasta, "fasta"):
            # Check if the ID in the header matches any in the list
            # We only need the second part of the header (e.g., "BnDHT1_076515" from ">g71571-T1 BnDHT1_076515")
            header_id = record.description.split()[1]
            if header_id in ids_to_extract:
                SeqIO.write(record, output_handle, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract sequences from a FASTA file based on a list of IDs.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-l", "--list", required=True, help="File containing list of IDs to extract")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file with extracted sequences")
    args = parser.parse_args()
    
    extract_fasta(args.input, args.list, args.output)
