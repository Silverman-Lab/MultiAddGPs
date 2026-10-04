# Run from the folder containing these three scripts.
# Read saved posterior draws and regenerate Web Figures 3-5. No model fitting.

library(ggplot2)
library(gridExtra)

project_dir <- "."
data_dir <- file.path(project_dir, "data/Data0_Simulation/data")
result_dir <- data_dir
output_dir <- file.path(project_dir, "figures")
source(file.path(project_dir, "simulation2_helper.R"))

# Data and saved results ------------------------------------------------

data_env <- new.env()
load(file.path(data_dir, "simulated_data2.RData"), envir = data_env)
sim <- data_env$simulated_data2
fit <- readRDS(file.path(result_dir, "Simulation_result2.rds"))
coordinates <- 1:3

# Linear batch-effect posterior draws.
linear_draws <- array(dim = c(sim$D - 1, sim$N, sim$samplesize))

for (s in seq_len(dim(linear_draws)[3])) {
  linear_draws[, , s] <- fit$mod$Lambda[[1]][, , s] %*% sim$X[1:2, ]
}

# In the saved simulation:
# f3 = periodic effect
# f2 = smooth trend
#
# In the fitted model:
# Lambda[[2]] = periodic effect
# Lambda[[3]] = smooth trend
components <- list(
  batch = list(
    NAddGP = center_over_time(fit$sample_posterior[[1]]),
    MultiAddGP = center_over_time(linear_draws),
    True = center_over_time(sim$f1X1)
  ),

  periodic = list(
    NAddGP = center_over_time(fit$sample_posterior[[2]]),
    MultiAddGP = center_over_time(fit$mod$Lambda[[2]]),
    True = center_over_time(sim$f3)
  ),

  trend = list(
    NAddGP = center_over_time(fit$sample_posterior[[3]]),
    MultiAddGP = center_over_time(fit$mod$Lambda[[3]]),
    True = center_over_time(sim$f2)
  ),

  total = list(
    NAddGP = fit$sample_posterior[[4]],
    MultiAddGP = fit$predicted,
    True = sim$f1X1 + sim$f2 + sim$f3
  )
)

# Figure styling --------------------------------------------------------

custom_colors <- c("#E43D40", "#0D698B", "#746C70", "#658EA9")

figure_theme <- theme_minimal() +
  theme(
    axis.title.y = element_text(size = 20, margin = margin(t = 15)),
    axis.title.x = element_text(size = 20, margin = margin(t = 15)),
    axis.text.x = element_text(angle = 0, hjust = 1, size = 17),
    axis.text.y = element_text(angle = 0, hjust = 1, size = 17),
    plot.title = element_text(size = 27, margin = margin(b = 20)),
    plot.subtitle = element_text(size = 20, margin = margin(b = 25))
  )

# Generate the three figures --------------------------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# Avoid creating an incidental Rplots.pdf while arranging the panels.
grDevices::pdf(NULL)

for (coordinate in coordinates) {
  summaries <- list(
    batch = summarize_component(components$batch, coordinate, sim$X1, "_f1"),
    periodic = summarize_component(components$periodic, coordinate, sim$X1, "_f2"),
    trend = summarize_component(components$trend, coordinate, sim$X1, "_f3"),
    total = summarize_component(components$total, coordinate, sim$X1)
  )

  # Panel A: batch effect.
  batch_plot <- component_plot(summaries$batch, "A",
    expression(F(X) == b[0] + b[1] * x^(Batch)), figure_theme, custom_colors, batch_effect = TRUE)

  # Panel B: periodic effect.
  periodic_plot <- component_plot(summaries$periodic, "B", expression(F(X) == f^(periodic)),
    figure_theme, custom_colors)

  # Panel C: smooth trend.
  trend_plot <- component_plot(summaries$trend, "C", expression(F(X) == f^(trend)),
    figure_theme, custom_colors)

  # Panel D: total fitted function.
  total_plot <- ggplot(summaries$total, aes(x = Date, y = mean)) +
    geom_ribbon(aes(ymin = p2.5, ymax = p97.5, fill = source), alpha = 0.2) +
    geom_ribbon(aes(ymin = p25, ymax = p75, fill = source), alpha = 0.5) +
    facet_wrap(~ batch, ncol = 1) +
    geom_line(aes(color = source, linetype = batch), linewidth = 0.7) +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = custom_colors) +
    figure_theme +
    theme(
      legend.title = element_blank(),
      strip.text = element_text(size = 20),
      legend.text = element_text(size = 17),
      legend.spacing.y = grid::unit(16, "pt"),
      legend.key.height = grid::unit(1.5, "cm")
    ) +
    labs(
      y = "F", x = "Hour", title = "D",
      subtitle = expression(F(X) == b[0] + b[1] * x^(Batch) + f^(periodic) + f^(trend))
    )

  # Preserve the original panel layout.
  layout_matrix <- rbind(c(1, 2, 4), c(3, 3, 4))

  figure <- arrangeGrob(batch_plot, trend_plot, periodic_plot, total_plot,
              layout_matrix = layout_matrix, widths = c(1.1, 1.3, 1.8))

  filename <- file.path(output_dir, sprintf("Web_Figure%d_corrected.png", coordinate + 2))
  ggsave(filename, plot = figure, width = 4893, height = 2477, units = "px", dpi = 300)
  ggsave(sub("[.]png$", ".pdf", filename), plot = figure,
         width = 4893 / 300, height = 2477 / 300)
  message("Saved ", filename)
}
grDevices::dev.off()

# RMSE and 95% coverage for this one dataset -----------------------------

metrics <- list()
for (component in names(components)) {
  for (method in c("NAddGP", "MultiAddGP")) {
    for (coordinate in coordinates) {
      posterior <- summarize_coordinate(components[[component]][[method]], coordinate, sim$X1)
      truth <- summarize_coordinate(components[[component]]$True, coordinate, sim$X1)$mean
      metrics[[length(metrics) + 1]] <- data.frame(
        component = component, method = method, coordinate = coordinate,
        RMSE = sqrt(mean((posterior$mean - truth)^2)),
        coverage = mean(truth >= posterior$p2.5 & truth <= posterior$p97.5)
      )
    }
  }
}
metrics <- do.call(rbind, metrics)
write.csv(metrics, file.path(result_dir, "single_dataset_metrics.csv"), row.names = FALSE)
