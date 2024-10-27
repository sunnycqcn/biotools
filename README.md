# biotools

python extract_gene_sequence.py -i fasta -l list -o extract.fa

sh run_blast.sh -d <database_file> -q <query_file> [-o <output_prefix>]

Rscript EDEGR_RNA_seq.R -i sample.list -n 3 -p 0.01 -f 0.05 -l 1 -c 10 -g comparisons.txt
Parameter Instructions
-i / --input:
Description: Path to the input count table or a file listing sample files (e.g., sample.list).
Required: Yes
Example: -i all_count_table.txt or -i sample.list
-n / --factor_num:
Description: The number of characters used to define groups from the sample names.
Default: 3
Example: -n 3 (assumes the first three characters of sample names define the group)
-p / --pvalue:
Description: P-value threshold for filtering the results of differential expression analysis.
Default: 0.05
Example: -p 0.01 (to filter results with a P-value less than or equal to 0.01)
-f / --fdr:
Description: False Discovery Rate (FDR) threshold for filtering the results.
Default: 0.05
Example: -f 0.05 (to filter results with FDR less than or equal to 0.05)
-l / --logFC:
Description: Log-fold change (logFC) threshold for filtering differentially expressed genes.
Default: 2
Example: -l 1 (to filter results with a logFC greater than or equal to 1)
-c / --cpm:
Description: Counts Per Million (CPM) threshold for filtering genes (CPM must be greater than or equal to this value in at least 2 samples).
Default: 10
Example: -c 10 (to keep genes with a CPM of at least 10)
-g / --compare_list:
Description: Path to a file listing specific comparisons to make. Each line should contain two groups (e.g., GroupA GroupB).
Default: NULL
Example: -g comparisons.txt (to specify particular groups for comparison)
