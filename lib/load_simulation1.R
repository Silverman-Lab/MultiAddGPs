# Load Simulation 1 dependencies and functions in dependency order.

run_dir <- Sys.getenv("MULTIADDGP_RUN_DIR", unset = "")
if (!nzchar(run_dir)) {
  stop("MULTIADDGP_RUN_DIR must be set before sourcing lib/load_simulation1.R")
}

source(file.path(run_dir, "lib", "loading_package.R"))
source(file.path(run_dir, "lib", "source.R"))
source(file.path(run_dir, "lib", "simulation_dgp.R"))
source(file.path(run_dir, "lib", "model_inputs.R"))
source(file.path(run_dir, "lib", "formatting.R"))
source(file.path(run_dir, "lib", "fit_multiaddgp.R"))
source(file.path(run_dir, "lib", "fit_lgpr.R"))
source(file.path(run_dir, "lib", "fit_gam.R"))
source(file.path(run_dir, "lib", "metrics.R"))
