# Aggregate per-seed metrics and runtimes.

args <- commandArgs(trailingOnly = TRUE)
run_dir <- normalizePath(
  if (length(args)) args[1] else getwd(),
  mustWork = TRUE
)

suppressPackageStartupMessages(library(dplyr))

models <- c("multiadd", "gam", "lgpr_gaussian", "lgpr_nb")

metric_rows <- list()
runtime_rows <- list()

for (model in models) {

  manifest <- read.csv(
    file.path(run_dir, "manifests", paste0(model, ".csv"))
  )

  paths <- file.path(
    run_dir,
    "results",
    model,
    manifest$cell_id,
    sprintf("result_seed%d.rds", manifest$seed)
  )

  if (!all(file.exists(paths))) {
    stop("Missing result: ", paths[which(!file.exists(paths))[1]])
  }

  for (i in seq_along(paths)) {

    result <- readRDS(paths[i])

    # Metrics
    m <- result$metrics
    m$model <- model
    m$D <- manifest$D[i]
    m$N <- manifest$N[i]
    m$seed <- manifest$seed[i]

    metric_rows[[length(metric_rows) + 1]] <- m

    # Runtime
    runtime <- result$runtime

    runtime_rows[[length(runtime_rows) + 1]] <- data.frame(
      model = model,
      D = manifest$D[i],
      N = manifest$N[i],
      seed = manifest$seed[i],
      observed_runtime_sec = runtime,
      estimated_total_runtime_sec =
        runtime * if (model == "multiadd") 1 else manifest$D[i]
    )
  }
}

metrics <- bind_rows(metric_rows)
runtimes <- bind_rows(runtime_rows)


# Summarize metrics.
metric_summary <- metrics |>
  group_by(model, D, N, component) |>
  summarise(
    n = n_distinct(seed),
    RMSE = mean(RMSE),
    Coverage = mean(Coverage),
    Avg_CI_Width = mean(Avg_CI_Width),
    .groups = "drop"
  )

metric_summary$cumulative_target <- ifelse(
  metric_summary$component == "Cumulative",
  "latent_CLR",
  NA
)

metric_summary$comparison_status <- "matched_target"

# NB cumulative CLR uses observed counts of the other taxa.
metric_summary$comparison_status[
  metric_summary$model == "lgpr_nb" &
    metric_summary$component == "Cumulative"
] <- "observed_other_taxa_plugin"


# Summarize runtime.
runtime_summary <- runtimes |>
  group_by(model, D, N) |>
  summarise(
    n = n_distinct(seed),
    mean_observed_runtime_sec = mean(observed_runtime_sec),
    mean_estimated_total_runtime_sec =
      mean(estimated_total_runtime_sec),
    .groups = "drop"
  )


# Save summaries.
out_dir <- file.path(run_dir, "table")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  metric_summary,
  file.path(out_dir, "metric_summary.csv"),
  row.names = FALSE
)

write.csv(
  runtime_summary,
  file.path(out_dir, "runtime_summary.csv"),
  row.names = FALSE
)

cat(
  "Saved",
  nrow(metric_summary), "metric rows and",
  nrow(runtime_summary), "runtime rows to",
  out_dir, "\n"
)
