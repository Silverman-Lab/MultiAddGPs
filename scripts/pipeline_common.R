Sys.setenv(MULTIADDGP_LIBRARY_ONLY = "1")

source(file.path(Sys.getenv("MULTIADDGP_RUN_DIR"), "scripts", "experiment_config.R"))
source(file.path(run_dir, "lib", "load_simulation1.R"))


# fun to read one row from a task csv
read_manifest_task <- function(manifest_path, task_id) {
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
  if (task_id < 1L || task_id > nrow(manifest)) {
    stop("task_id ", task_id, " is outside manifest range 1..", nrow(manifest))
  }
  manifest[task_id, , drop = FALSE]
}

# construct a filename for a simualted dataset
dataset_path <- function(D, N, seed) {
  file.path(data_dir, sprintf("D%d_N%d", D, N), sprintf("sim_seed%d.rds", seed))
}

# construct a filename for a model result
model_result_path <- function(model, D, N, seed) {
  file.path(result_dir, model, sprintf("D%d_N%d", D, N),sprintf("result_seed%d.rds", seed))
}

# # Save safely: write to a temporary file first, then move it into place.
atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp.", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) {
    stop("Could not atomically move result into place: ", path)
  }
}

# Check that the RDS file exists, is readable, and contains the expected objects.
is_valid_rds <- function(path, required_names = character()) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    return(FALSE)
  }
  object <- tryCatch(readRDS(path), error = function(e) NULL)
  !is.null(object) && all(required_names %in% names(object))
}

# Record software versions so the simulation can be reproduced later.
software_versions <- function() {
  packages <- c(
    "lgpr", "rstan", "mgcv", "dplyr", "ParBayesianOptimization",
    "fido"
  )
  versions <- vapply(packages, function(package) {
    if (requireNamespace(package, quietly = TRUE)) {
      as.character(utils::packageVersion(package))
    } else {
      NA_character_
    }
  }, character(1))
  list(R = R.version.string, packages = versions)
}
