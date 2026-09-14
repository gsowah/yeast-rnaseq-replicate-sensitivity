# Running the full-data DESeq2 analysis comparing knockout with wild type.

suppressPackageStartupMessages({
  library(DESeq2)
})

metadata_file <- "data/metadata/schurch_metadata.csv"
counts_file <- "data/processed/schurch_counts.csv"
output_dir <- "results/full_data"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.csv(
  metadata_file,
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

count_matrix <- read.csv(
  counts_file,
  row.names = 1,
  check.names = FALSE
)

count_matrix <- as.matrix(count_matrix)

# If statement for checking count matrix and metadata contain the same samples
if (!setequal(colnames(count_matrix), rownames(metadata))) {
  stop("Sample names in the count matrix and metadata don't match.")
}

# Putting metadata rows in the same order as count-matrix columns
metadata <- metadata[colnames(count_matrix), , drop = FALSE]

# wild type as the reference condition
metadata$condition <- factor(
  metadata$condition,
  levels = c("wildtype", "knockout")
)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) > 0, ]

# differential-expression analysis with DESEq2
dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "knockout", "wildtype"),
  alpha = 0.05
)

res <- res[order(res$padj, na.last = TRUE), ]

results_table <- cbind(
  gene_id = rownames(res),
  as.data.frame(res)
)

write.csv(
  results_table,
  file.path(output_dir, "deseq2_knockout_vs_wildtype.csv"),
  row.names = FALSE
)

normalized_counts <- counts(dds, normalized = TRUE)

write.csv(
  data.frame(
    gene_id = rownames(normalized_counts),
    normalized_counts
  ),
  file.path(output_dir, "normalized_counts.csv"),
  row.names = FALSE
)
