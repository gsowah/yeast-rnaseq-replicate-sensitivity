# This script performs repeated random subsampling across different replicate
# numbers, re-runs DESeq2 for each randomly selected set of samples.
# 3 to 10 replicates per condition, 100 runs per replicate number.
# Measure recovery of significant genes from the full-data analysis

suppressPackageStartupMessages({
    library(DESeq2)
})

# Main parameters

replicate_numbers <- 3:10
runs_per_n <- 100
base_seed <- 123
significance_threshold <- 0.05


metadata_file <- "data/metadata/schurch_metadata.csv"
counts_file <- "data/processed/schurch_counts.csv"
full_results_file <- file.path(
    "results",
    "full_data",
    "deseq2_knockout_vs_wildtype.csv"
)


output_dir <- file.path(
    "results",
    "subsampling",
    "repeated_runs"
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

# Validation of analysis parameters
# QC ensuring parameters are valid before 800 analysis run

if (
    length(replicate_numbers) == 0 ||
    anyNA(replicate_numbers) ||
    any(replicate_numbers < 2) ||
    any(replicate_numbers != round(replicate_numbers))
) {
    stop("Replicate numbers must be whole numbers of at least 2.")
}

if (
    length(runs_per_n) != 1 ||
    is.na(runs_per_n) ||
    runs_per_n < 1 ||
    runs_per_n != round(runs_per_n)
) {
    stop("runs_per_n must be one positive whole number.")
}

if (
    length(significance_threshold) != 1 ||
    is.na(significance_threshold) ||
    significance_threshold <= 0 ||
    significance_threshold >= 1
) {
    stop("The significance threshold must be between 0 and 1.")
}

# Wildtype as reference

metadata$condition <- factor(
    metadata$condition,
    levels = c("wildtype", "knockout")
)

wildtype_samples <- rownames(metadata)[
    metadata$condition == "wildtype"
]

knockout_samples <- rownames(metadata)[
    metadata$condition == "knockout"
]

if (
    max(replicate_numbers) > length(wildtype_samples) ||
    max(replicate_numbers) > length(knockout_samples)
) {
    stop("Not enough samples for the requested replicate numbers.")
}

# retaining genes less than 0.05 pertaining to adjusted p values

full_significant_genes <- full_results$gene_id[
    !is.na(full_results$padj) &
        full_results$padj < significance_threshold
]

# storing the repeated-run results
# 800 reruns since 8 (total rep numbers) * 100 (number of runs)

total_runs <- length(replicate_numbers) * runs_per_n # replicate numbers from 3-10

run_summaries <- vector(
    mode = "list",
    length = total_runs
)

result_index <- 1

# for loops for repeated re-runs
# Run repeated subsampling analysis

for (n in replicate_numbers) {

    cat("Starting analyses for n = ", n, "\n", sep = "")

    for (run_number in seq_len(runs_per_n)) {

        run_seed <- base_seed + result_index - 1
        set.seed(run_seed)

        selected_wildtype <- sample(
            wildtype_samples,
            size = n,
            replace = FALSE
        )

        selected_knockout <- sample(
            knockout_samples,
            size = n,
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

        # subsampled DESeq2 model

        subsample_dds <- DESeqDataSetFromMatrix(
            countData = selected_counts,
            colData = selected_metadata,
            design = ~ condition
        )

        # removing genes with 0 reads

        subsample_dds <- subsample_dds[
            rowSums(counts(subsample_dds)) > 0,
        ]

        subsample_dds <- DESeq(
            subsample_dds,
            quiet = TRUE
        )

        subsample_results <- results(
            subsample_dds,
            contrast = c(
                "condition",
                "knockout",
                "wildtype"
            ),
            alpha = significance_threshold
        )

        # identifying and comparing the significant genes

        subsample_significant_genes <- rownames(subsample_results)[
            !is.na(subsample_results$padj) &
                subsample_results$padj < significance_threshold
        ]

        recovered_genes <- intersect(
            full_significant_genes,
            subsample_significant_genes
        )

        subsample_only_genes <- setdiff(
            subsample_significant_genes,
            full_significant_genes
        )

        recovery_percentage <- (
            length(recovered_genes) /
                length(full_significant_genes)
        ) * 100

        # Run summary

        run_summaries[[result_index]] <- data.frame(
            replicate_number = n,
            run_number = run_number,
            random_seed = run_seed,
            total_selected_samples = length(selected_samples),
            full_significant_genes = length(full_significant_genes),
            subsample_significant_genes =
                length(subsample_significant_genes),
            recovered_genes = length(recovered_genes),
            subsample_only_genes = length(subsample_only_genes),
            recovery_percentage = recovery_percentage
        )

        result_index <- result_index + 1
    }
}

# Adding results from all the runs

all_runs_table <- do.call(
    rbind,
    run_summaries
)

rownames(all_runs_table) <- NULL

# Making a summary table to add all results

summary_table <- do.call(
    rbind,
    lapply(
        split(
            all_runs_table,
            all_runs_table$replicate_number
        ),
        function(run_data) {
            data.frame(
                replicate_number =
                    unique(run_data$replicate_number),
                runs = nrow(run_data),
                mean_subsample_significant_genes =
                    mean(run_data$subsample_significant_genes),
                mean_recovered_genes =
                    mean(run_data$recovered_genes),
                mean_subsample_only_genes =
                    mean(run_data$subsample_only_genes),
                mean_recovery_percentage =
                    mean(run_data$recovery_percentage),
                standard_deviation =
                    sd(run_data$recovery_percentage),
                median_recovery_percentage =
                    median(run_data$recovery_percentage),
                minimum_recovery_percentage =
                    min(run_data$recovery_percentage),
                maximum_recovery_percentage =
                    max(run_data$recovery_percentage)
            )
        }
    )
)

rownames(summary_table) <- NULL


write.csv(
    all_runs_table,
    file.path(output_dir, "all_runs.csv"), # containing 800 individualruns
    row.names = FALSE
)

write.csv(
    summary_table,
    file.path(output_dir, "summary_by_replicate_number.csv"), # containing the 8 summary rows
    row.names = FALSE
)