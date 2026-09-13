# Script will randomly select three wild-type and three knockout samples
# Rerun DESeq2 on the six selected samples
# Compare gene recovery with the full 86-sample analysis

suppressPackageStartupMessages({
    library(DESeq2)
})

# Main parameters
random_seed <- 123
replicates_per_condition <- 3
significance_threshold <- 0.05

metadata_file <- "data/metadata/schurch_metadata.csv"
counts_file <- "data/processed/schurch_counts.csv"
full_results_file <- paste0(
    "results/full_data/",
    "deseq2_knockout_vs_wildtype.csv"
)

output_dir <- file.path(
    "results",
    "subsampling",
    paste0("n", replicates_per_condition, "_single_run")
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


metadata <- read.csv(
    metadata_file,
    row.names = 1,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

counts <- read.csv(
    counts_file,
    row.names = 1,
    check.names = FALSE
)

counts <- as.matrix(counts)

full_results <- read.csv(
    full_results_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# Random selection of samples per condition

set.seed(random_seed)

wildtype_samples <- rownames(metadata)[
    metadata$condition == "wildtype"
]

knockout_samples <- rownames(metadata)[
    metadata$condition == "knockout"
]

selected_wildtype <- sample(
    wildtype_samples,
    size = replicates_per_condition,
    replace = FALSE
)

selected_knockout <- sample(
    knockout_samples,
    size = replicates_per_condition,
    replace = FALSE
)

selected_samples <- c(
    selected_wildtype,
    selected_knockout
)

selected_metadata <- metadata[
    selected_samples,
    ,
    drop = FALSE
]

selected_counts <- counts[
    ,
    selected_samples,
    drop = FALSE
]

# run *subsampled DESeq2 model

subsample_dds <- DESeqDataSetFromMatrix(
    countData = selected_counts,
    colData = selected_metadata,
    design = ~ condition
)

# 0 read genes removal

subsample_dds <- subsample_dds[
    rowSums(counts(subsample_dds)) > 0,
]

subsample_dds <- DESeq(subsample_dds)


subsample_results <- results(
    subsample_dds,
    contrast = c("condition", "knockout", "wildtype"),
    alpha = significance_threshold
)

subsample_results <- subsample_results[
    order(subsample_results$padj, na.last = TRUE),
]

# Identify significant genes using adjusted p values

full_significant_genes <- full_results$gene_id[
    !is.na(full_results$padj) &
        full_results$padj < significance_threshold
]

subsample_significant_genes <- rownames(subsample_results)[
    !is.na(subsample_results$padj) &
        subsample_results$padj < significance_threshold
]

# ComparING significant-gene sets

recovered_genes <- intersect(
    full_significant_genes,
    subsample_significant_genes
)

subsample_only_genes <- setdiff(
    subsample_significant_genes,
    full_significant_genes
)

if (length(full_significant_genes) == 0) {
    stop("The full-data analysis contains no significant genes.")
}

recovery_percentage <- (
    length(recovered_genes) /
        length(full_significant_genes)
) * 100

# output tables

subsample_results_table <- as.data.frame(subsample_results)
subsample_results_table$gene_id <- rownames(subsample_results_table)

subsample_results_table <- subsample_results_table[
    ,
    c(
        "gene_id",
        setdiff(colnames(subsample_results_table), "gene_id")
    )
]

selected_samples_table <- data.frame(
    sample_id = selected_samples,
    condition = as.character(selected_metadata$condition),
    stringsAsFactors = FALSE
)

recovered_genes_table <- data.frame(
    gene_id = recovered_genes,
    stringsAsFactors = FALSE
)

subsample_only_genes_table <- data.frame(
    gene_id = subsample_only_genes,
    stringsAsFactors = FALSE
)

recovery_summary_table <- data.frame(
    random_seed = random_seed,
    replicates_per_condition = replicates_per_condition,
    total_selected_samples = length(selected_samples),
    full_significant_genes = length(full_significant_genes),
    subsample_significant_genes = length(subsample_significant_genes),
    recovered_genes = length(recovered_genes),
    subsample_only_genes = length(subsample_only_genes),
    recovery_percentage = recovery_percentage
)


write.csv(
    selected_samples_table,
    file.path(output_dir, "selected_samples.csv"),
    row.names = FALSE
)

write.csv(
    subsample_results_table,
    file.path(output_dir, "deseq2_knockout_vs_wildtype.csv"),
    row.names = FALSE
)

write.csv(
    recovered_genes_table,
    file.path(output_dir, "recovered_genes.csv"),
    row.names = FALSE
)

write.csv(
    subsample_only_genes_table,
    file.path(output_dir, "subsample_only_genes.csv"),
    row.names = FALSE
)

write.csv(
    recovery_summary_table,
    file.path(
        output_dir,
        "single_run_recovery_summary.csv"
    ),
    row.names = FALSE
)
