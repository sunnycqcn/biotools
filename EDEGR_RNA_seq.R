#!/usr/bin/env Rscript

# Function to check and install required packages
check_and_install_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg, repos = "http://cran.us.r-project.org")
      library(pkg, character.only = TRUE)
    }
  }
}

# List of necessary packages
packages <- c("edgeR", "optparse", "pheatmap", "RColorBrewer", "ggplot2")
check_and_install_packages(packages)

# Define command-line arguments
library(optparse)
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
    
    count_data[[basename(sample_file)]] <- sample_counts  
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
data <- as.matrix(data)

# Filter genes based on CPM threshold (CPM >= opt$cpm in at least 2 samples)
library(edgeR)
d <- DGEList(counts = data)
d <- calcNormFactors(d)
cps <- cpm(d)
keep <- rowSums(cps >= cpm_thresh) >= 2
data <- data[keep, ]
grp <- as.factor(substr(colnames(data), 1, factor_num))

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
      
      contrast_name <- paste0(group1, "-vs-", group2)
      con <- makeContrasts(contrasts = paste0("grp", group1, "-grp", group2), levels = mm)
      
      cat("Running differential expression for contrast:", contrast_name, "\n")
      
      f <- glmFit(d, mm)
      lrt <- glmLRT(f, contrast = con)
      
      tt <- topTags(lrt, n = Inf)$table
      filtered_tt <- tt[tt$PValue <= p_value_thresh & tt$FDR <= fdr_thresh & abs(tt$logFC) >= logFC_thresh, ]
      
      if (nrow(filtered_tt) > 0) {
        filtered_tt$gene_id <- rownames(filtered_tt)
        colnames(filtered_tt)[colnames(filtered_tt) == "logFC"] <- paste0("logFC_", contrast_name)
      }
      
      combined_results[[contrast_name]] <- filtered_tt
      all_gene_ids <- unique(c(all_gene_ids, filtered_tt$gene_id))
      
      write.table(filtered_tt, file = paste0(contrast_name, "_filtered.xls"), sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }
  
  combined_table <- data.frame(gene_id = all_gene_ids)
  for (contrast_name in names(combined_results)) {
    contrast_data <- combined_results[[contrast_name]]
    if (nrow(contrast_data) > 0) {
      combined_table <- merge(combined_table, contrast_data[, c("gene_id", paste0("logFC_", contrast_name))], by = "gene_id", all.x = TRUE)
    }
  }
  
  write.table(combined_table, file = "combined_logFC_results.xls", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

# Run DGE
run_dge(d, mm, levels(grp), p_value_thresh, fdr_thresh, logFC_thresh, compare_list_file)

# Plot pairwise sample comparison
plot_pairwise_comparison <- function(d) {
  plotMDS(d, main = "Pairwise Sample Comparison (MDS Plot)")
  dev.copy(png, "pairwise_sample_comparison.png")
  dev.off()
}

# BCV Plot
plot_bcv <- function(d) {
  plotBCV(d, main = "Biological Coefficient of Variation (BCV)")
  dev.copy(png, "bcv_plot.png")
  dev.off()
}

# Mean-Variance plot
plot_mean_variance <- function(d) {
  plotMeanVar(d, main = "Mean-Variance Plot")
  dev.copy(png, "mean_variance_plot.png")
  dev.off()
}

# PCA Plot
plot_pca <- function(d) {
  pca <- prcomp(t(cpm(d)))
  df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], group = d$samples$group)
  
  ggplot(df, aes(PC1, PC2, color = group)) +
    geom_point(size = 3) +
    labs(title = "PCA Plot", x = "PC1", y = "PC2") +
    theme_minimal()
  
  ggsave("pca_plot.png")
}

# Heatmap using group means for DEGs
plot_heatmap_group_means <- function(d, degs) {
  group_means <- sapply(levels(d$samples$group), function(g) rowMeans(d$counts[degs, d$samples$group == g, drop = FALSE]))
  log_group_means <- log2(group_means + 1)
  
  png("heatmap_DEGs_group_means.png")
  pheatmap(log_group_means, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = TRUE,
           color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))
  dev.off()
}

# Get DEGs for heatmap
degs <- rownames(d$counts)[rowSums(cpm(d) > 1) >= 2]
plot_heatmap_group_means(d, degs)

# Call plots
plot_pairwise_comparison(d)
plot_bcv(d)
plot_mean_variance(d)
plot_pca(d)

cat("All analyses and plots have been generated.\n")

