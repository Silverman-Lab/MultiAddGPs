args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript generate_dataset_task.R MANIFEST.csv TASK_ID")
}

run_dir <- normalizePath(getwd(), mustWork = TRUE)
Sys.setenv(MULTIADDGP_RUN_DIR = run_dir)
source(file.path(run_dir, "scripts", "pipeline_common.R"))

task <- read_manifest_task(args[[1]], as.integer(args[[2]]))
path <- dataset_path(task$D, task$N, task$seed)

if (is_valid_rds(path, c("Y", "X", "D", "N"))) {
  cat("SKIP valid dataset:", path, "\n")
  quit(save = "no")
}

dataset <- sim(
  D = task$D,
  N = task$N,
  hyper_params = hyper_params,
  multi.depth = multi_depth,
  seed = task$seed
)

dataset$experiment_metadata <- list(
  cell_id = task$cell_id,
  seed_index = task$seed_index,
  seed = task$seed
)

atomic_save_rds(dataset, path)
cat("Saved dataset:", path, "\n")
