# Plot metric and runtime summaries.

args <- commandArgs(trailingOnly = TRUE)
run_dir <- normalizePath(
  if (length(args)) args[1] else getwd(),
  mustWork = TRUE
)

extra_libs <- Sys.getenv("MULTIADDGP_R_LIBS", "")
if (nzchar(extra_libs)) {
  .libPaths(c(strsplit(extra_libs, .Platform$path.sep)[[1]], .libPaths()))
}

library(ggplot2)
library(dplyr)
library(patchwork)

summary_dir <- file.path(run_dir, "table")
out_dir <- file.path(run_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# Plot settings.
model_labels <- c(
  multiadd = "MultiAddGPs",
  lgpr_gaussian = "Lgpr_Gaussian",
  lgpr_nb = "Lgpr_nb",
  gam = "GAM"
)

model_levels <- unname(model_labels)

model_colors <- c(
  MultiAddGPs = "#0a27ca",
  Lgpr_Gaussian = "#e69f00",
  Lgpr_nb = "#e65c00",
  GAM = "#009E73"
)

components <- c(
  "Batch", "General", "Intervention", "Cumulative"
)

format_models <- function(data) {
  data$model[data$model == "lgpr_nb_observed_other_taxa"] <- "lgpr_nb"

  data %>%
    mutate(
      model = factor(
        recode(model, !!!model_labels),
        levels = model_levels
      ),
      component = if ("component" %in% names(data)) {
        factor(component, levels = components)
      } else {
        NULL
      }
    )
}

base_theme <- theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )


# Metric plots.
metrics <- read.csv(
  file.path(summary_dir, "metric_summary.csv")
) %>%
  format_models()

create_metric_grid <- function(data, x_var, x_label, title) {
  data$x_axis <- factor(data[[x_var]])

  rmse_plot <- ggplot(
    data,
    aes(x = x_axis, y = RMSE, color = model, group = model)
  ) +
    guides(color = "none") +
    facet_grid(. ~ component, scales = "free_y") +
    geom_line(linewidth = 0.45, na.rm = TRUE) +
    geom_point(size = 1.35, na.rm = TRUE) +
    scale_color_manual(values = model_colors, drop = FALSE) +
    labs(title = title, x = NULL, y = "RMSE") +
    base_theme +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      legend.position = "none"
    )

  coverage_plot <- ggplot(
    data,
    aes(x = x_axis, y = Coverage, fill = model)
  ) +
    facet_grid(. ~ component) +
    geom_col(position = "dodge") +
    geom_hline(
      yintercept = 0.95,
      linetype = "dashed",
      alpha = 0.5
    ) +
    scale_fill_manual(values = model_colors, drop = FALSE) +
    coord_cartesian(ylim = c(0, 1.1)) +
    labs(x = x_label, y = "Coverage", fill = "Model") +
    base_theme

  rmse_plot / coverage_plot
}

n_plot <- create_metric_grid(
  filter(metrics, D == 3),
  "N",
  "Sample size (N)",
  "(A). D fixed, N varying"
) &
  theme(legend.position = "none")

d_plot <- create_metric_grid(
  filter(metrics, N == 50),
  "D",
  "Taxa count (D)",
  "(B). N fixed, D varying"
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

final_metrics_plot <- (
  wrap_elements(full = n_plot) /
    wrap_elements(full = d_plot)
) +
  plot_layout(heights = c(1, 1))

ggsave(
  file.path(out_dir, "final_metrics_plot.pdf"),
  final_metrics_plot,
  width = 8.5,
  height = 10
)

ggsave(
  file.path(out_dir, "final_metrics_plot.png"),
  final_metrics_plot,
  width = 8.5,
  height = 10,
  dpi = 300
)


# Runtime plots.
runtime <- read.csv(
  file.path(summary_dir, "runtime_summary.csv")
) %>%
  mutate(
    model = factor(
      recode(model, !!!model_labels),
      levels = model_levels
    )
  )

create_runtime_plot <- function(data, x_var, x_label, y_label = NULL) {
  data$x_axis <- factor(data[[x_var]])

  ggplot(
    data,
    aes(
      x = x_axis,
      y = mean_estimated_total_runtime_sec,
      color = model,
      group = model
    )
  ) +
    geom_line(linewidth = 0.45, na.rm = TRUE) +
    geom_point(size = 1.6, na.rm = TRUE) +
    scale_y_log10(
      labels = scales::label_number(accuracy = 1)
    ) +
    scale_color_manual(
      values = model_colors,
      name = "Model",
      drop = FALSE
    ) +
    labs(x = x_label, y = y_label) +
    theme_bw()
}

p_run_d <- create_runtime_plot(
  filter(runtime, N == 50),
  "D",
  "Taxa (D)",
  "Total Runtime (s)"
)

p_run_n <- create_runtime_plot(
  filter(runtime, D == 3),
  "N",
  "Samples (N)"
)

final_runtime_plot <- (p_run_d | p_run_n) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  file.path(out_dir, "final_runtime_plot.pdf"),
  final_runtime_plot,
  width = 7,
  height = 4
)

ggsave(
  file.path(out_dir, "final_runtime_plot.png"),
  final_runtime_plot,
  width = 7,
  height = 4,
  dpi = 300
)

cat("Saved final plots to", out_dir, "\n")
