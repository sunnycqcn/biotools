###python rename_gff.py -i input.gff -o output.gff

import importlib
import subprocess
import sys
import re
import argparse

def install(module_name):
    """Install a module using pip."""
    subprocess.check_call([sys.executable, "-m", "pip", "install", module_name])

def check_and_install_modules():
    """Check for required modules and install them if not found."""
    required_modules = []
    
    # Check for the argparse module
    try:
        import argparse
    except ImportError:
        required_modules.append("argparse")
    
    # Add other modules here if necessary
    # e.g., 
    # try:
    #     import some_module
    # except ImportError:
    #     required_modules.append("some_module")
    
    # Install any missing modules
    for module in required_modules:
        print(f"Installing module: {module}...")
        install(module)

def rename_gff(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        gene_ids = {}  # Store mappings from old to new gene IDs
        mRNA_ids = {}  # Store mappings for mRNA and associated IDs
        parent_ids = {}  # Store mappings for parent IDs
        
        chromosome_counter = {}  # To count the number of genes per chromosome

        # First pass: collect mappings and rename gene IDs
        for line in infile:
            # Skip header lines
            if line.startswith("##"):
                outfile.write(line)
                continue
            
            parts = line.strip().split("\t")
            if len(parts) < 9:
                outfile.write(line)
                continue
            
            # Extract chromosome number
            chromosome = parts[0]
            prefix = ''.join(filter(str.isalpha, chromosome))  # Get the string prefix
            chromosome_number = re.sub(r'\D', '', chromosome)  # Extract digits
            
            # Initialize counter for new chromosome
            if chromosome not in chromosome_counter:
                chromosome_counter[chromosome] = 1
            
            # Rename gene IDs
            original_id_match = re.search(r'ID=(SH_\d+)', parts[8])
            if original_id_match:
                original_id = original_id_match.group(1)
                # Format new gene ID
                new_gene_id = f"{prefix}{int(chromosome_number):02d}G{chromosome_counter[chromosome]:06d}"
                gene_ids[original_id] = new_gene_id  # Store the mapping
                
                # Update the gene line with the new ID
                parts[8] = re.sub(r'ID=[^;]+', f'ID={new_gene_id}', parts[8])
                outfile.write("\t".join(parts) + "\n")
                chromosome_counter[chromosome] += 1  # Increment gene count for this chromosome
                continue
            
            # Rename mRNA IDs and update their parents
            mRNA_id_match = re.search(r'ID=(g\d+-T\d+)', parts[8])
            if mRNA_id_match:
                mRNA_id = mRNA_id_match.group(1)
                # Get the new mRNA ID based on the corresponding gene ID
                if original_id in gene_ids:
                    new_mRNA_id = f"{gene_ids[original_id]}_{mRNA_id.split('-')[-1]}"
                    mRNA_ids[mRNA_id] = new_mRNA_id  # Store the mapping for later reference
                    parts[8] = re.sub(r'ID=[^;]+', f'ID={new_mRNA_id}', parts[8])
                    parent_id_match = re.search(r'Parent=(g\d+-T\d+)', parts[8])
                    if parent_id_match:
                        parent_id = parent_id_match.group(1)
                        if parent_id in mRNA_ids:
                            parts[8] = re.sub(r'Parent=[^;]+', f'Parent={mRNA_ids[parent_id]}', parts[8])
                outfile.write("\t".join(parts) + "\n")
                continue
            
            # Rename child features (exons, UTRs, CDS) based on mRNA mapping
            if 'Parent=' in parts[8]:
                parent_id_match = re.search(r'Parent=(g\d+-T\d+)', parts[8])
                if parent_id_match:
                    parent_id = parent_id_match.group(1)
                    if parent_id in mRNA_ids:
                        parts[8] = re.sub(r'Parent=[^;]+', f'Parent={mRNA_ids[parent_id]}', parts[8])
                        
                        # Change IDs based on mRNA ID format
                        new_feature_id = f"{mRNA_ids[parent_id]}{re.search(r'\.(.*)', parts[8]).group(0)}"
                        parts[8] = re.sub(r'ID=[^;]+', f'ID={new_feature_id}', parts[8])
            
            # Write the modified line to the output file
            outfile.write("\t".join(parts) + "\n")

if __name__ == "__main__":
    check_and_install_modules()  # Check and install required modules

    # Argument parser for input and output file handling
    parser = argparse.ArgumentParser(description="Rename GFF file IDs based on specific rules.")
    parser.add_argument('-i', '--input', required=True, help='Input GFF file')
    parser.add_argument('-o', '--output', required=True, help='Output GFF file')
    args = parser.parse_args()
    
    rename_gff(args.input, args.output)
