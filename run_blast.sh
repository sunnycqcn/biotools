# Default parameters
num_threads=20
evalue="1e-200"
max_targets=5

# Parse command-line arguments
while getopts "d:q:o:" opt; do
    case $opt in
        d) db_file="$OPTARG" ;;        # Database file
        q) query_file="$OPTARG" ;;     # Query file
        o) output_prefix="$OPTARG" ;;  # Output prefix (optional)
        *) echo "Invalid option"; exit 1 ;;
    esac
done

# Check for required arguments
if [[ -z "$db_file" || -z "$query_file" ]]; then
    echo "Usage: sh blast.sh -d <database_file> -q <query_file> [-o <output_prefix>]"
    exit 1
fi

# Function to detect sequence type (nucl or prot)
detect_sequence_type() {
    local file=$1
    # Check if the sequence contains any protein-specific residues
    if grep -iq "[EFILPQZ]" "$file"; then
        echo "prot"
    else
        echo "nucl"
    fi
}

# Detect sequence types for query and database
db_type=$(detect_sequence_type "$db_file")
query_type=$(detect_sequence_type "$query_file")

# Determine the correct BLAST tool based on detected types
if [[ "$db_type" == "nucl" && "$query_type" == "nucl" ]]; then
    blast_tool="blastn"
elif [[ "$db_type" == "prot" && "$query_type" == "prot" ]]; then
    blast_tool="blastp"
elif [[ "$db_type" == "nucl" && "$query_type" == "prot" ]]; then
    blast_tool="tblastn"
elif [[ "$db_type" == "prot" && "$query_type" == "nucl" ]]; then
    blast_tool="blastx"
else
    echo "Error: Could not determine the appropriate BLAST tool"
    exit 1
fi

echo "Detected database type: $db_type"
echo "Detected query type: $query_type"
echo "Using BLAST tool: $blast_tool"

# Define output file prefix
output_prefix="${output_prefix:-blast_output}"

# Create the BLAST database if needed
makeblastdb -in "$db_file" -dbtype "$db_type" -parse_seqids -out db_blast

# Run the selected BLAST tool with specified formats
declare -A formats=( [6]="6" [7]="7" [10]="10" [0]="0" )
for fmt in "${!formats[@]}"; do
    output_file="${output_prefix}.${formats[$fmt]}.out"
    $blast_tool -query "$query_file" -db db_blast -evalue "$evalue" -max_target_seqs "$max_targets" -outfmt "$fmt" -out "$output_file" -num_threads "$num_threads"
done

echo "BLAST search completed with $blast_tool. Results saved to ${output_prefix}.[6|7|10|0].out"

