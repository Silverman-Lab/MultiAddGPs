run_dir <- normalizePath(
  Sys.getenv("MULTIADDGP_RUN_DIR", getwd()),
  mustWork = TRUE
)

# Define paths
config_dir     <- file.path(run_dir, "config")
manifest_dir   <- file.path(run_dir, "manifests")
data_dir       <- file.path(run_dir, "data")
result_dir     <- file.path(run_dir, "results")


# Simulation grid setting 
simulation_grid <- data.frame(
  D = c(rep(3L, 7), 10L, 50L, 100L, 300L, 500L),
  N = c(20L, 50L, 70L, 100L, 300L, 500L, 1000L, rep(50L, 5))
)
simulation_grid$cell_id <- sprintf(
  "D%d_N%d", simulation_grid$D, simulation_grid$N
)
simulation_grid <- simulation_grid[, c("cell_id", "D", "N")]


num_seeds <- 100L
multi_depth <- 2000L
hyper_params <- c(2, 3, 1, 0.5)

# detemine D and N for each model
model_is_enabled <- function(model, D, N) {
  switch(
    model,
    multiadd = TRUE,
    gam = TRUE,
    lgpr_gaussian =
      (D == 3L && N %in% c(20L, 50L, 70L, 100L)) ||
      (N == 50L && D == 10L),
    lgpr_nb =
      D == 3L && N %in% c(20L, 50L),
    stop("Unknown model: ", model)
  )
}

# set the # of seeds for each setting, currently is all the same
model_seed_count <- function(model, D, N) num_seeds

# This is only for lgpr, specify the mcmc setting for lgpr_gassuian and lgpr_nb
# as nb required run substantially larger iteration.
sampler_config <- function(model, N) {
  switch(
    model,
    lgpr_gaussian = list(iter = 6500L, chains = 4L, cores = 4L),
    lgpr_nb = list(
      iter = if (N == 20L) 6500L else 8000L,
      chains = 4L,
      cores = 4L
    ),
    list(iter = NA_integer_, chains = NA_integer_, cores = 1L)
  )
}

# set the flag thresold for lgpr later
ess_threshold <- 2000
rhat_threshold <- 1.01
