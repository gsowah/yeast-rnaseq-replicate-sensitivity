# Making plots from DESeq full data analysis 

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

dds_file <- "results/full_data/full_data_dds.rds"
figure_dir <- "results/figures/full_data_qc"
output_dir <- "results/full_data"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE
)

# Load the fitted DESeq2 object
dds <- readRDS(dds_file)

# Transform counts for visualisation plots like pca and heatmap
vsd <- vst(dds, blind = FALSE)

condition_colours <- c(
  wildtype = "#0072B2",
  knockout = "#D55E00"
)

# PCA plot
pca_data <- plotPCA(
  vsd,
  intgroup = "condition",
  returnData = TRUE
)

percent_variance <- round(
  100 * attr(pca_data, "percentVar"),
  digits = 1
)

pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    colour = condition
  )
) +
  geom_point(size = 3, alpha = 0.8) +
  scale_colour_manual(values = condition_colours) +
  labs(
    title = "PCA of full Schurch RNA-seq dataset",
    x = paste0("PC1: ", percent_variance[1], "% variance"),
    y = paste0("PC2: ", percent_variance[2], "% variance"),
    colour = "Condition"
  ) +
  theme_classic(base_size = 12)

ggsave(
  filename = file.path(figure_dir, "pca_condition.png"),
  plot = pca_plot,
  width = 7,
  height = 5,
  units = "in",
  dpi = 300
)

# Sample-distance heatmap
sample_distances <- dist(t(assay(vsd)))
sample_distance_matrix <- as.matrix(sample_distances)

rownames(sample_distance_matrix) <- colnames(vsd)
colnames(sample_distance_matrix) <- colnames(vsd)

sample_annotations <- as.data.frame(
  colData(vsd)[, "condition", drop = FALSE]
)

annotation_colours <- list(
  condition = condition_colours
)

png(
  filename = file.path(figure_dir, "sample_distance_heatmap.png"),
  width = 2400,
  height = 2200,
  res = 300
)

pheatmap(
  sample_distance_matrix,
  annotation_col = sample_annotations,
  annotation_row = sample_annotations,
  annotation_colors = annotation_colours,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Sample distances after variance-stabilising transformation"
)

dev.off()

# Library-size plot
library_sizes <- colSums(counts(dds))

library_data <- data.frame(
  sample = names(library_sizes),
  library_size = as.numeric(library_sizes),
  condition = colData(dds)$condition
)

library_plot <- ggplot(
  library_data,
  aes(
    x = reorder(sample, library_size),
    y = library_size,
    fill = condition
  )
) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = condition_colours) +
  scale_y_continuous(
    labels = scales::comma
  ) +
  labs(
    title = "Raw library sizes",
    x = "Sample",
    y = "Total read count",
    fill = "Condition"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 5)
  )

ggsave(
  filename = file.path(figure_dir, "library_sizes.png"),
  plot = library_plot,
  width = 8,
  height = 11,
  units = "in",
  dpi = 300
)

# MA plot
png(
  filename = file.path(figure_dir, "ma_plot.png"),
  width = 2100,
  height = 1800,
  res = 300
)

plotMA(
  dds,
  alpha = 0.05,
  ylim = c(-10, 10),
  main = "Knockout versus wild-type"
)

dev.off()

# Dispersion plot
png(
  filename = file.path(figure_dir, "dispersion_plot.png"),
  width = 2100,
  height = 1800,
  res = 300
)

plotDispEsts(
  dds,
  main = "DESeq2 dispersion estimates"
)

dev.off()

pca_output <- data.frame(
  sample = rownames(pca_data),
  pca_data,
  row.names = NULL
)

write.csv(
  pca_output,
  file.path(output_dir, "pca_coordinates.csv"),
  row.names = FALSE
)

cat("\nPCA and quality-control outputs created successfully.\n")
cat("Figure directory:", figure_dir, "\n")
cat("PCA coordinates saved to:", output_dir, "\n")
cat(
  "PC1 and PC2 variance:",
  paste0(percent_variance, "%", collapse = ", "),
  "\n"
)
