# MultiAddGP hyperparameter sensitivity: prepare, fit, and plot RMSE and coverage.
# The model implementation remains in lib/fit_multiaddgp.R.

parameter_names <- c("sigma_treat", "rho_treat", "sigma_general", "rho_general")

component_levels <- c("Batch", "General", "Intervention", "Cumulative")

setting_levels <- c("Baseline", as.vector(t(outer(parameter_names, c("-10%", "+10%"), paste))))


# Setup ------------------------------------------------------------

initialize_sensitivity <- function(root) {
  Sys.setenv(MULTIADDGP_RUN_DIR = normalizePath(root, mustWork = TRUE))
  source(file.path(Sys.getenv("MULTIADDGP_RUN_DIR"), "scripts", "pipeline_common.R"))
}


select_baselines <- function(D, N, seed_count = 10L) {
  manifest <- read.csv(file.path(manifest_dir, "multiadd.csv"), stringsAsFactors = FALSE)

  manifest <- manifest[
    manifest$D == D & manifest$N == N,
    ,
    drop = FALSE
  ]

  manifest <- manifest[
    order(manifest$seed_index),
    ,
    drop = FALSE
  ]

  head(manifest, seed_count)
}


# Paths ------------------------------------------------------------

sensitivity_tag <- function(D, N, seed_count) {
  sprintf("D%d_N%d_%dseed", D, N, seed_count)
}


sensitivity_manifest_path <- function(D, N, seed_count) {
  file.path(
    manifest_dir,
    paste0("multiadd_sensitivity_", sensitivity_tag(D, N, seed_count), ".csv")
  )
}


sensitivity_result_path <- function(task) {
  file.path(
    result_dir,
    "multiadd_sensitivity",
    task$cell_id,
    sprintf("result_seed%d_%s_%s.rds", task$seed, task$parameter, task$perturbation)
  )
}


# Prepare sensitivity tasks ----------------------------------------

prepare_sensitivity <- function(D, N, seed_count = 10) {
  source_manifest <- select_baselines(D, N, seed_count)

  perturbations <- data.frame(
    perturbation = c("minus10", "plus10"),
    factor = c(0.9, 1.1),
    stringsAsFactors = FALSE
  )

  rows <- list()

  for (i in seq_len(nrow(source_manifest))) {
    task <- source_manifest[i, , drop = FALSE]

    for (param_index in seq_along(parameter_names)) {
      for (j in seq_len(nrow(perturbations))) {

        rows[[length(rows) + 1L]] <- data.frame(
          model = "multiadd",
          cell_id = task$cell_id,
          D = task$D,
          N = task$N,
          seed_index = task$seed_index,
          seed = task$seed,
          param_index = param_index,
          parameter = parameter_names[[param_index]],
          perturbation = perturbations$perturbation[[j]],
          factor = perturbations$factor[[j]],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (!length(rows)) {
    stop("No baseline tasks selected.")
  }

  manifest <- do.call(rbind, rows)
  manifest$task_id <- seq_len(nrow(manifest))

  manifest <- manifest[, c(
    "task_id", "model", "cell_id", "D", "N",
    "seed_index", "seed", "param_index",
    "parameter", "perturbation", "factor"
  )]

  out <- sensitivity_manifest_path(D, N, seed_count)

  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  write.csv(manifest, out, row.names = FALSE)

  cat("Wrote", nrow(manifest), "sensitivity tasks to", out, "\n")

  invisible(out)
}


# Fit one sensitivity task -----------------------------------------

fit_sensitivity_task <- function(manifest_path, task_id) {
  task <- read_manifest_task(manifest_path, task_id)

  dataset <- dataset_path(task$D, task$N, task$seed)
  baseline_path <- model_result_path("multiadd", task$D, task$N, task$seed)

  baseline <- readRDS(baseline_path)

  base_parameters <- baseline$diagnostics$bo_parameters
  names(base_parameters) <- parameter_names

  perturbed_parameters <- base_parameters
  perturbed_parameters[[task$param_index]] <-
    perturbed_parameters[[task$param_index]] * task$factor

  output <- sensitivity_result_path(task)

  if (is_valid_rds(output, c("metadata", "estimates", "metrics", "runtime"))) {
    cat("SKIP valid sensitivity result:", output, "\n")
    return(invisible(output))
  }

  sim <- readRDS(dataset)
  truth <- extract_truth(sim)

  fit <- fit_multiaddgps(sim, perturbed_parameters, task$seed)

  metrics <- evaluate_methods(setNames(list(fit$estimates), "multiadd_sensitivity"), truth)

  result <- list(
    metadata = list(
      model = "multiadd_sensitivity",
      baseline_model = "multiadd",
      cell_id = task$cell_id,
      D = task$D,
      N = task$N,
      seed_index = task$seed_index,
      seed = task$seed,
      parameter = task$parameter,
      perturbation = task$perturbation,
      factor = task$factor,
      completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    ),
    baseline_parameters = base_parameters,
    perturbed_parameters = perturbed_parameters,
    estimates = fit$estimates,
    metrics = metrics,
    runtime = fit$runtime,
    software = software_versions()
  )

  atomic_save_rds(result, output)

  cat("Saved sensitivity result:", output, "\n")

  invisible(output)
}


# Collect completed results ----------------------------------------

collect_sensitivity <- function(manifest, baseline_manifest) {
  records <- list()

  for (i in seq_len(nrow(baseline_manifest))) {
    task <- baseline_manifest[i, , drop = FALSE]

    path <- model_result_path("multiadd", task$D, task$N, task$seed)

    object <- readRDS(path)

    parameters <- object$diagnostics$bo_parameters
    names(parameters) <- parameter_names

    records[[length(records) + 1L]] <- list(
      object = object,
      task = task,
      setting = "Baseline",
      parameter = "baseline",
      perturbation = "baseline",
      factor = 1,
      parameters = parameters
    )
  }

  for (i in seq_len(nrow(manifest))) {
    task <- manifest[i, , drop = FALSE]
    path <- sensitivity_result_path(task)

    object <- readRDS(path)

    records[[length(records) + 1L]] <- list(
      object = object,
      task = task,
      setting = paste(task$parameter, ifelse(task$perturbation == "minus10", "-10%", "+10%")),
      parameter = task$parameter,
      perturbation = task$perturbation,
      factor = task$factor,
      parameters = object$perturbed_parameters
    )
  }

  if (!length(records)) {
    stop("No completed results available.")
  }

  records
}


# Summarize multiple seeds -----------------------------------------

summarize_seeds <- function(records) {
  metric_rows <- lapply(records, function(record) {
    metrics <- record$object$metrics

    metrics$model <- "multiadd"
    metrics$setting <- record$setting
    metrics$parameter <- record$parameter
    metrics$perturbation <- record$perturbation
    metrics$factor <- record$factor
    metrics$D <- record$task$D
    metrics$N <- record$task$N
    metrics$seed_index <- record$task$seed_index
    metrics$seed <- record$task$seed

    metrics
  })

  metrics <- do.call(rbind, metric_rows)

  metrics$setting <- factor(metrics$setting, levels = setting_levels)

  metrics$component <- factor(metrics$component, levels = component_levels)


  summarize_metric <- function(values) {
    n <- sum(is.finite(values))
    mean_value <- mean(values, na.rm = TRUE)
    sd_value <- sd(values, na.rm = TRUE)
    se_value <- sd_value / sqrt(n)

    data.frame(
      n = n,
      mean = mean_value,
      sd = sd_value,
      se = se_value,
      ci_lower = mean_value - 1.96 * se_value,
      ci_upper = mean_value + 1.96 * se_value
    )
  }


  metric_summary <- do.call(
    rbind,
    lapply(
      split(metrics, interaction(metrics$setting, metrics$component, drop = TRUE)),
      function(z) {
        rmse <- summarize_metric(z$RMSE)
        coverage <- summarize_metric(z$Coverage)

        data.frame(
          setting = z$setting[1],
          parameter = z$parameter[1],
          perturbation = z$perturbation[1],
          factor = z$factor[1],
          component = z$component[1],

          n_seeds = rmse$n,

          RMSE = rmse$mean,
          RMSE_sd = rmse$sd,
          RMSE_se = rmse$se,
          RMSE_ci_lower = rmse$ci_lower,
          RMSE_ci_upper = rmse$ci_upper,

          Coverage = coverage$mean,
          Coverage_sd = coverage$sd,
          Coverage_se = coverage$se,
          Coverage_ci_lower = pmax(0, coverage$ci_lower),
          Coverage_ci_upper = pmin(1, coverage$ci_upper),

          stringsAsFactors = FALSE
        )
      }
    )
  )

  metric_summary$setting <- factor(metric_summary$setting, levels = setting_levels)

  metric_summary$component <- factor(metric_summary$component, levels = component_levels)

  metric_summary$is_baseline <-
    metric_summary$setting == "Baseline"

  cat("Completed metric rows:", nrow(metrics), "\n")
  cat("Completed settings:", length(unique(metrics$setting)), "of", length(setting_levels), "\n")

  metric_summary
}


# Plot multiple-seed summary ---------------------------------------

plot_seed_summary <- function(metric_summary, out_dir, D, N, seed_count) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  metric_summary$setting <- factor(metric_summary$setting, levels = setting_levels)

  metric_summary$component <- factor(metric_summary$component, levels = component_levels)

  baseline_summary <- metric_summary[
    metric_summary$is_baseline,
    ,
    drop = FALSE
  ]

  theme_sensitivity <- theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )


  p_rmse <- ggplot(metric_summary, aes(y = setting)) +
    geom_vline(data = baseline_summary, aes(xintercept = RMSE), linetype = "dotted", color = "black", linewidth = 0.45) +
    geom_segment(aes(x = RMSE_ci_lower, xend = RMSE_ci_upper, yend = setting), color = "grey35", linewidth = 0.65) +
    geom_point(aes(x = RMSE, shape = is_baseline), color = "black", fill = "white", size = 2.3) +
    scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 16)) +
    facet_wrap(~ component, scales = "free_x", ncol = 2) +
    labs(title = sprintf("MultiAddGP hyperparameter sensitivity (D = %d, N = %d, %d seeds)", D, N, seed_count), x = "RMSE", y = "Perturbation setting") +
    theme_sensitivity

  ggsave(
    file.path(out_dir, "multiadd_sensitivity_rmse_ci.png"),
    p_rmse,
    width = 10,
    height = 6,
    dpi = 300
  )

  ggsave(file.path(out_dir, "multiadd_sensitivity_rmse_ci.pdf"), p_rmse, width = 10, height = 6)


  p_coverage <- ggplot(metric_summary, aes(y = setting)) +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = "grey45") +
    geom_vline(data = baseline_summary, aes(xintercept = Coverage), linetype = "dotted", color = "black", linewidth = 0.45) +
    geom_segment(aes(x = Coverage_ci_lower, xend = Coverage_ci_upper, yend = setting), color = "grey35", linewidth = 0.65) +
    geom_point(aes(x = Coverage, shape = is_baseline), color = "black", fill = "white", size = 2.3) +
    scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 16)) +
    facet_wrap(~ component, ncol = 2) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(title = sprintf("MultiAddGP hyperparameter sensitivity coverage (D = %d, N = %d, %d seeds)", D, N, seed_count), x = "Coverage", y = "Perturbation setting") +
    theme_sensitivity

  ggsave(
    file.path(out_dir, "multiadd_sensitivity_coverage_ci.png"),
    p_coverage,
    width = 10,
    height = 6,
    dpi = 300
  )

  ggsave(
    file.path(out_dir, "multiadd_sensitivity_coverage_ci.pdf"),
    p_coverage,
    width = 10,
    height = 6
  )

  cat("Saved sensitivity plots to", out_dir, "\n")
}


# Report -----------------------------------------------------------

report_sensitivity <- function(D, N, seed_count = 10L) {
  manifest <- read.csv(sensitivity_manifest_path(D, N, seed_count), stringsAsFactors = FALSE)

  baseline_manifest <- select_baselines(D, N, seed_count)

  records <- collect_sensitivity(manifest, baseline_manifest)

  out_dir <- file.path(run_dir, "figures", "sensitivity", sensitivity_tag(D, N, seed_count))

  metric_summary <- summarize_seeds(records)

  plot_seed_summary(metric_summary, out_dir, D, N, seed_count)

  invisible(out_dir)
}


# Command line interface -------------------------------------------

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) < 2L) {
    stop(
      "Usage: multiadd_sensitivity.R MODE RUN_DIR ",
      "[D N SEED_COUNT]; ",
      "fit MODE: RUN_DIR MANIFEST TASK_ID"
    )
  }

  mode <- match.arg(args[[1]], c("prepare", "fit", "report", "all"))

  initialize_sensitivity(args[[2]])

  if (mode == "fit") {
    if (length(args) != 4L) {
      stop("fit requires RUN_DIR MANIFEST TASK_ID")
    }

    return(invisible(fit_sensitivity_task(args[[3]], as.integer(args[[4]]))))
  }

  D <- if (length(args) >= 3L) {
    as.integer(args[[3]])
  } else {
    3L
  }

  N <- if (length(args) >= 4L) {
    as.integer(args[[4]])
  } else {
    50L
  }

  seed_count <- if (length(args) >= 5L) {
    as.integer(args[[5]])
  } else {
    10L
  }

  if (anyNA(c(D, N, seed_count)) || any(c(D, N, seed_count) < 1L)) {
    stop("D, N, and SEED_COUNT must be positive integers")
  }

  if (mode %in% c("prepare", "all")) {
    path <- prepare_sensitivity(D, N, seed_count)
  }

  if (mode == "all") {
    manifest <- read.csv(path)

    for (task_id in seq_len(nrow(manifest))) {
      fit_sensitivity_task(path, task_id)
    }
  }

  if (mode %in% c("report", "all")) {
    report_sensitivity(D, N, seed_count)
  }

  invisible(NULL)
}


if (sys.nframe() == 0L) {
  main()
}
