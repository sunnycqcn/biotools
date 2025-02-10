#!/usr/bin/env Rscript

# Load necessary libraries
check_and_install_packages <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if (length(new_packages)) {
    options(repos = c(CRAN = "https://cloud.r-project.org/"))  # Set CRAN mirror
    install.packages(new_packages)
  }
}

# List of required packages
required_packages <- c("edgeR", "optparse", "reshape2", "ggplot2", "pheatmap")
check_and_install_packages(required_packages)

# Load libraries
library(edgeR)
library(optparse)
library(reshape2)
library(ggplot2)
library(pheatmap)

# Define command-line arguments
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "Path to input count table (e.g., all_count_table.txt or sample.list)", metavar = "character"),
  make_option(c("-n", "--factor_num"), type = "integer", default = 3, 
              help = "Number of characters used to define groups (default: 3)"),
  make_option(c("-p", "--pvalue"), type = "numeric", default = 0.05, 
              help = "P-value threshold for filtering results (default: <= 0.05)"),
  make_option(c("-f", "--fdr"), type = "numeric", default = 0.05, 
              help = "FDR threshold for filtering results (default: <= 0.05)"),
  make_option(c("-l", "--logFC"), type = "numeric", default = 2, 
              help = "logFC threshold for filtering results (default: >= 2)"),
  make_option(c("-c", "--cpm"), type = "numeric", default = 10, 
              help = "CPM threshold for filtering genes (default: CPM >= 10 in at least 2 samples)"),
  make_option(c("-g", "--compare_list"), type = "character", default = NULL, 
              help = "Path to comparison list file (format: group1 group2 per line)", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if input file is provided
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file is required. Use -i to specify the input file.", call. = FALSE)
}

# Read input arguments
input_file <- opt$input
factor_num <- opt$factor_num
p_value_thresh <- opt$pvalue
fdr_thresh <- opt$fdr
logFC_thresh <- opt$logFC
cpm_thresh <- opt$cpm
compare_list_file <- opt$compare_list

# Function to generate all_count_table.txt from sample.list
generate_count_table <- function(sample_list_file, output_file) {
  sample_files <- readLines(sample_list_file)  # Read the list of sample files
  count_data <- list()  # Initialize a list to store count data
  
  for (sample_file in sample_files) {
    sample_counts <- read.table(sample_file, header = FALSE, stringsAsFactors = FALSE)
    colnames(sample_counts) <- c("gene_id", basename(sample_file))
    
    sample_counts <- sample_counts[!sample_counts$gene_id %in% 
      c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique"), ]
    
    sample_counts[[2]] <- as.numeric(gsub("[^0-9]", "0", sample_counts[[2]]))
    
    if (nrow(sample_counts) == 0) {
      cat("No counts remaining after filtering for sample:", sample_file, "\n")
      next
    }
    
    count_data[[basename(sample_file)]] <- sample_counts  
  }
  
  if (length(count_data) == 0) {
    stop("No valid sample files provided. Please check your sample list.", call. = FALSE)
  }
  
  all_counts <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), count_data)
  
  write.table(all_counts, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
}

if (grepl("sample.list$", input_file)) {
  all_count_table_file <- "all_count_table.txt"
  generate_count_table(input_file, all_count_table_file)
  input_file <- all_count_table_file
}

data <- read.table(input_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE)

if (is.null(data) || nrow(data) == 0) {
  stop("Error: The input count data is empty. Please check your input data.", call. = FALSE)
}

cat("Input count data preview:\n")
print(head(data))
cat("Data dimensions:", dim(data), "\n")

data <- as.matrix(data)

if (!is.numeric(data)) {
  data[!sapply(data, is.numeric)] <- 0
  data <- apply(data, 2, function(x) as.numeric(as.character(x)))
}

if (!is.matrix(data) || length(dim(data)) < 2) {
  stop("Error: The data is not a matrix or does not have at least two dimensions.", call. = FALSE)
}

if (any(rowSums(data) == 0)) {
  cat("Removing genes with all zero counts.\n")
  data <- data[rowSums(data) > 0, ]
}

if (nrow(data) == 0) {
  stop("Error: The count data is empty after filtering. Please check your input data.", call. = FALSE)
}

# Filter genes based on CPM threshold (CPM >= opt$cpm in at least 2 samples)
d <- DGEList(counts = data)
d <- calcNormFactors(d)
cps <- cpm(d)

keep <- rowSums(cps >= cpm_thresh) >= 2
data <- data[keep, ]

if (nrow(data) == 0) {
  stop("Error: The count data is empty after filtering for CPM. Please check your input data.", call. = FALSE)
}

grp <- as.factor(substr(colnames(data), 1, factor_num))
cat("Group summary:\n")
print(table(grp))

d <- DGEList(counts = data, group = grp)
d <- calcNormFactors(d)

mm <- model.matrix(~ -1 + grp)
d <- estimateGLMCommonDisp(d, mm)
d <- estimateGLMTrendedDisp(d, mm)
d <- estimateGLMTagwiseDisp(d, mm)

# Function to run DGE and combine results
run_dge <- function(d, mm, grp_levels, p_value_thresh, fdr_thresh, logFC_thresh, compare_list = NULL) {
  combined_results <- list()
  all_gene_ids <- character(0)
  
  # If comparison list is provided, read it
  if (!is.null(compare_list)) {
    compare_pairs <- read.table(compare_list, header = FALSE, stringsAsFactors = FALSE)
  }
  
  for (i in 1:(length(grp_levels) - 1)) {
    for (j in (i + 1):length(grp_levels)) {
      group1 <- grp_levels[i]
      group2 <- grp_levels[j]
      
      # If comparison list is provided, skip groups not in list
      if (!is.null(compare_list) && !any(compare_pairs$V1 == group1 & compare_pairs$V2 == group2)) {
        next
      }
      
      contrast_name <- paste0(group1, "-vs-", group2)
      con <- makeContrasts(contrasts = paste0("grp", group1, "-grp", group2), levels = mm)
      
      cat("Running differential expression for contrast:", contrast_name, "\n")
      
      f <- glmFit(d, mm)
      lrt <- glmLRT(f, contrast = con)
      
      tt <- topTags(lrt, n = Inf)$table
      
      # Rename columns for clarity
      colnames(tt) <- paste0(colnames(tt), "_", contrast_name)

      filtered_tt <- tt[tt$PValue <= p_value_thresh & tt$FDR <= fdr_thresh & abs(tt$logFC) >= logFC_thresh, ]
      
      if (nrow(filtered_tt) > 0) {
        filtered_tt$gene_id <- rownames(filtered_tt)
      } else {
        filtered_tt <- data.frame(gene_id = character(0), logFC = numeric(0), PValue = numeric(0), FDR = numeric(0))
      }
      
      combined_results[[contrast_name]] <- filtered_tt
      all_gene_ids <- unique(c(all_gene_ids, filtered_tt$gene_id))
    }
  }
  
  combined_table <- data.frame(gene_id = all_gene_ids)
  for (contrast in names(combined_results)) {
    combined_table <- merge(combined_table, combined_results[[contrast]], by = "gene_id", all.x = TRUE)
  }
  
  write.table(combined_table, file = "combined_logFC_results.xls", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

# Run the DGE analysis
run_dge(d, mm, levels(grp), p_value_thresh, fdr_thresh, logFC_thresh, compare_list_file)

# Function for pairwise comparison plot
pairwise_comparison_plot <- function(d) {
  comparison <- as.data.frame(d$counts)
  comparison <- comparison[rowSums(comparison) > 0, ]
  
  comparison_long <- reshape2::melt(comparison, variable.name = "Sample", value.name = "Count")
  
  # Ensure Gene column exists for facetting
  comparison_long$Gene <- rownames(comparison)
  
  ggplot(comparison_long, aes(x = Sample, y = Count)) +
    geom_boxplot() +
    facet_wrap(~ Gene) +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("Pairwise Sample Comparison Plot") +
    ylab("Count") +
    xlab("Sample")
}

# Create pairwise comparison plot
pairwise_comparison_plot(d)

# Function for BCV plot
bcv_plot <- function(d) {
  plotBCV(d)
  title("Biological Coefficient of Variation Plot")
}

# Create BCV plot
bcv_plot(d)

# Function for mean-variance plot
mean_variance_plot <- function(d) {
  plotMeanVar(d, show.text = TRUE, pch = 20, col = ifelse(d$tagwise.dispersion < 0.5, "blue", "red"))
  title("Mean-Variance Plot")
}

# Create mean-variance plot
mean_variance_plot(d)

# Function for PCA plot
pca_plot <- function(d) {
  pca <- prcomp(t(d$counts), center = TRUE, scale. = TRUE)
  pca_df <- as.data.frame(pca$x)
  pca_df$group <- d$samples$group
  
  ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 3) +
    ggtitle("PCA Plot") +
    xlab("Principal Component 1") +
    ylab("Principal Component 2")
}

# Create PCA plot
pca_plot(d)

# Function to plot heatmap of DEGs
heatmap_deg <- function(d, p_value_thresh, fdr_thresh) {
  results <- read.table("combined_logFC_results.xls", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  degs <- results[results$FDR <= fdr_thresh & results$PValue <= p_value_thresh, ]
  
  if (nrow(degs) == 0) {
    stop("No DEGs found with the specified thresholds.", call. = FALSE)
  }
  
  heatmap_data <- d$counts[rownames(d$counts) %in% degs$gene_id, ]
  pheatmap(heatmap_data, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, 
           show_colnames = TRUE, main = "Heatmap of DEGs")
}

# Create heatmap of DEGs
heatmap_deg(d, p_value_thresh, fdr_thresh)

cat("Differential expression analysis completed.\n")
