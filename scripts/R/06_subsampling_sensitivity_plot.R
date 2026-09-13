# relationship between replicate number and gene recovery
# Uses the summary produced by the repeated subsampling analysis

suppressPackageStartupMessages({
    library(ggplot2)
})


summary_file <- file.path(
    "results",
    "subsampling",
    "repeated_runs",
    "summary_by_replicate_number.csv"
)

output_file <- file.path(
    "results",
    "figures",
    "subsampling",
    "sensitivity_vs_replicate_number.png"
)

dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
)


summary_table <- read.csv(
    summary_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# plot

sensitivity_plot <- ggplot(
    summary_table,
    aes(
        x = replicate_number,
        y = mean_recovery_percentage
    )
) +
    geom_errorbar(
        aes(
            ymin = pmax(
                0,
                mean_recovery_percentage - standard_deviation
            ),
            ymax = pmin(
                100,
                mean_recovery_percentage + standard_deviation
            )
        ),
        width = 0.15,
        linewidth = 0.7,
        colour = "#0072B2"
    ) +
    geom_line(
        linewidth = 1,
        colour = "#0072B2"
    ) +
    geom_point(
        size = 3,
        colour = "#0072B2"
    ) +
    scale_x_continuous(
        breaks = summary_table$replicate_number
    ) +
    scale_y_continuous(
        limits = c(0, 100),
        breaks = seq(0, 100, by = 10)
    ) +
    labs(
        title = "Gene recovery increases with replicate number",
        subtitle = "Mean ± standard deviation across 100 random subsamples",
        x = "Replicates per condition",
        y = "Recovery of full-data significant genes (%)"
    ) +
    theme_minimal(
        base_size = 12
    ) +
    theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(
            face = "bold"
        )
    )

ggsave(
    filename = output_file,
    plot = sensitivity_plot,
    width = 8,
    height = 5,
    units = "in",
    dpi = 300,
    bg = "white"
)