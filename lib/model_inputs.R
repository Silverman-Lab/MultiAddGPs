# Model input construction for Simulation 1.
#
# These helpers convert a simulated dataset into the data frames used by GAM
# and LGPR model fits.

# Normalize dataset columns to the raw-time design expected by MultiAddGP.
normalize_simulation_design <- function(sim) {
  if ("treatment" %in% rownames(sim$X)) {
    stopifnot(!("time" %in% rownames(sim$X)),
              is.numeric(sim$X1$time), length(sim$X1$time) == ncol(sim$X))
    sim$X["treatment", ] <- sim$X1$time
    rownames(sim$X)[rownames(sim$X) == "treatment"] <- "time"
  }
  if ("treatment" %in% names(sim$X1)) sim$X1$treatment <- NULL
  sim
}

make_gam_lgpr_nb_data <- function(sim, taxon_idx = 2) {
  stopifnot(is.list(sim), !is.null(sim$X), !is.null(sim$Y), !is.null(sim$X1))
  stopifnot(taxon_idx >= 1, taxon_idx <= nrow(sim$Y))

  X <- sim$X
  N <- ncol(X)
  X1 <- as.data.frame(t(X))
  # Construct the GAM/LGPR input data.
  data <- cbind(X1, Y = sim$Y[taxon_idx, ]) %>%
    mutate(Batch = as.factor(ifelse(batchBatch2 == 1, "Batch1", "Batch2")),
           ID = factor("001"),
           time = seq(0, N - 1, by = 1),
           Date2 = Date) %>%
    dplyr::select(-batchBatch2)
}

make_gam_lgpr_Gaussian_data <- function(sim, taxon_idx = 2) {
  stopifnot(is.list(sim), !is.null(sim$X), !is.null(sim$Y), !is.null(sim$X1))
  stopifnot(taxon_idx >= 1, taxon_idx <= nrow(sim$Y))

  X <- sim$X
  N <- ncol(X)
  X1 <- as.data.frame(t(X))
  # Construct the GAM/LGPR input data.
  data <- cbind(X1,  clr_Y = clr_array(sim$Y+0.65, parts = 1)[taxon_idx, ]) %>%
    mutate(Batch = as.factor(ifelse(batchBatch2 == 1, "Batch1", "Batch2")),
           ID = factor("001"),
           time = seq(0, N - 1, by = 1),
           Date2 = Date) %>%
    dplyr::select(-batchBatch2)
}
