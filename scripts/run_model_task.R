# Run one model task from the manifest.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript run_model_task.R MANIFEST.csv TASK_ID")
}

run_dir <- normalizePath(getwd(), mustWork = TRUE)
Sys.setenv(MULTIADDGP_RUN_DIR = run_dir)
source(file.path(run_dir, "scripts", "pipeline_common.R"))

# Load task and simulated dataset.
task <- read_manifest_task(args[1], as.integer(args[2]))
dataset <- dataset_path(task$D, task$N, task$seed)

if (!is_valid_rds(dataset, c("Y", "X", "D", "N"))) {
  stop("Missing or invalid dataset: ", dataset)
}

# Set output path.
output <- model_result_path(
  task$model, task$D, task$N, task$seed
)

suffix <- Sys.getenv("MODEL_RESULT_SUFFIX", "")

if (nzchar(suffix)) {output <- sub("\\.rds$", paste0("_", suffix, ".rds"), output)}

if (is_valid_rds(output,c("metadata", "estimates", "metrics", "runtime"))) {
  cat("SKIP valid model result:", output, "\n")
  quit(save = "no")
}

# Load data and truth.
set.seed(task$seed)
sim <- readRDS(dataset)
truth <- extract_truth(sim)

is_lgpr <- task$model %in% c("lgpr_gaussian", "lgpr_nb")

# Fit model.
if (task$model == "multiadd") {
  fit <- fit_multiadd_pipeline(sim, task$seed)

} else if (task$model == "gam") {
  data <- make_gam_lgpr_Gaussian_data(sim, taxon_idx = 2)
  fit <- fit_gam(data, seq_len(task$N))
  fit$diagnostics <- list()

} else if (is_lgpr) {
  fit <- fit_lgpr_pipeline(
    sim,
    task$model,
    task$iter,
    task$chains,
    task$cores,
    task$seed
  )

} else {
  stop("Unsupported model: ", task$model)
}

# Evaluate model.
metrics <- evaluate_methods(
  setNames(list(fit$estimates), task$model),
  truth
)

flagged <- is_lgpr && lgpr_is_flagged(fit$diagnostics)

sampler_control <- if (is.null(fit$sampler_control)) {
  list()
} else {
  fit$sampler_control
}

# Save results.
result <- list(
  metadata = list(
    model = task$model,
    cell_id = task$cell_id,
    D = task$D,
    N = task$N,
    seed_index = task$seed_index,
    seed = task$seed,
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    diagnostic_flag = flagged,
    cumulative_target = "latent_CLR"
  ),
  estimates = fit$estimates,
  summary_spec = fit$summary_spec,
  metrics = metrics,
  runtime = fit$runtime,
  diagnostics = fit$diagnostics,
  sampler = c(
    list(
      iter = task$iter,
      chains = task$chains,
      cores = task$cores
    ),
    sampler_control
  ),
  software = software_versions()
)

atomic_save_rds(result, output)

# Save Stan fit when LGPR diagnostics are flagged.
if (is_lgpr && flagged) {
  stanfit_path <- sub("\\.rds$", "_stanfit.rds", output)
  atomic_save_rds(fit$fit, stanfit_path)
}

cat("Saved model result:", output, "\n")
cat("Diagnostic flag:", flagged, "\n")