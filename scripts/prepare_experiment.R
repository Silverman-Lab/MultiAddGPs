# Load experiment configuration.
args <- commandArgs(trailingOnly = TRUE)
run_dir <- normalizePath(
  if (length(args)) args[1] else getwd(),
  mustWork = TRUE
)

Sys.setenv(MULTIADDGP_RUN_DIR = run_dir)
source(file.path(run_dir, "scripts", "experiment_config.R"))

# Create directories.
dirs <- c(config_dir, manifest_dir, data_dir, result_dir)
lapply(dirs,dir.create,recursive = TRUE,showWarnings = FALSE)


# Load fixed seeds.
seeds <- read.csv(file.path(config_dir, "seeds.csv"))
stopifnot(nrow(seeds) == num_seeds)
seeds$seed_index <- as.integer(seeds$seed_index)
seeds$seed <- as.integer(seeds$seed)


# Use the simulation grid defined in experiment_config.R.
grid <- simulation_grid

# create the grid cells with seeds together
cross_manifest <- function(grid) {
  out <- merge(grid, seeds, by = NULL)
  out[order(out$D, out$N, out$seed_index), ]
}

# Dataset manifest.
dataset_manifest <- cross_manifest(grid)
dataset_manifest$task_id <- seq_len(nrow(dataset_manifest))
dataset_manifest <- dataset_manifest[
  , c("task_id", "cell_id", "D", "N", "seed_index", "seed")
]

# save the dataset manifest in a csv file
write.csv(
  dataset_manifest,
  file.path(manifest_dir, "datasets.csv"),
  row.names = FALSE
)

# Model-specific manifests.
make_model_manifest <- function(model) {
  enabled <- mapply(model_is_enabled, model,grid$D,grid$N)
  out <- cross_manifest(grid[enabled, , drop = FALSE])
  # Keep only the required number of seeds for each model/grid cell.
  seed_limits <- mapply(
    model_seed_count,
    model,
    out$D,
    out$N
  )
  out <- out[out$seed_index <= seed_limits, ]
  out$model <- model
  # Sampler settings.
  configs <- lapply(out$N, function(N) sampler_config(model, N))

  out$iter <- vapply(configs, `[[`, integer(1), "iter")
  out$chains <- vapply(configs, `[[`, integer(1), "chains")
  out$cores <- vapply(configs, `[[`, integer(1), "cores")

  out$task_id <- seq_len(nrow(out))

  out[, c(
    "task_id", "model", "cell_id", "D", "N",
    "seed_index", "seed", "iter", "chains", "cores"
  )]
}

# Write one manifest for each model with specified D,N and MCMC setting
models <- c("multiadd","gam","lgpr_gaussian","lgpr_nb")
for (model in models) {
  write.csv(
    make_model_manifest(model),
    file.path(manifest_dir, paste0(model, ".csv")),
    row.names = FALSE
  )
}

# Summary.
cat(
  "Prepared experiment in:", run_dir,
  "\nGrid cells:", nrow(grid),
  "\nDatasets:", nrow(dataset_manifest),
  "\nFixed seeds:", nrow(seeds), "\n"
)
