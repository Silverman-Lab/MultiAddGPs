# MultiAddGP fitting for Simulation 1.
#
# These functions fit the MultiAddGP model and compute the marginal-likelihood
# score used during hyperparameter search.

fit_multiaddgps_for_mml <- function(sim, hyper_params, seed,
                                   samplesize = 300,
                                   max_iter = 200,
                                   jitter = 1e-4) {

  # Normalize the design matrix before fitting.
  sim <- normalize_simulation_design(sim)

  Gamma_general <- function(X) SE(X, hyper_params[3], hyper_params[4], X_row = 3, jitter = 1e-6)
  Gamma_intervention <- function(X) intervention_kernel(
    X, sigma = hyper_params[1], rho = hyper_params[2], X_row = "time", jitter = 1e-6
  )

  linear_component_theta <- matrix(0, sim$D - 1, 2)
  linear_component_gamma <- diag(2) + 1e-6

  Theta_kernel <- function(X) matrix(0, sim$D - 1, ncol(X))
  Theta <- list(linear_component_theta, Theta_kernel, Theta_kernel)
  Gamma <- list(linear_component_gamma, Gamma_general, Gamma_intervention)

  fit <- fido::basset(
    sim$Y, sim$X,
    upsilon = sim$D + 100,
    Theta = Theta,
    Gamma = Gamma,
    linear = c(1:2),
    n_samples = samplesize,
    verbose = FALSE,
    seed = seed + 100,
    jitter = jitter,
    max_iter = max_iter,
    init = alr_array(sim$Y + 0.65, parts = 1)
  )

  fit
}


MLL_multiaddgps <- function(sim,
                            sigma_treat, rho_treat,
                            sigma_general, rho_general,
                            seed = 893,
                            use_priors = TRUE,
                            prior_lambda = 1,
                            fast_samplesize = 300,
                            fast_max_iter = 200) {

  hyper_params <- c(sigma_treat, rho_treat, sigma_general, rho_general)

  fit <- fit_multiaddgps_for_mml(
    sim = sim,
    hyper_params = hyper_params,
    seed = seed,
    samplesize = fast_samplesize,
    max_iter = fast_max_iter
  )

  logML <- fit$logMarginalLikelihood

  if (!use_priors) {
    return(list(Score = logML))
  }

  # --- Example priors (tune to match your paper) ---
  # keep these weak at first; prior_lambda=1 is usually plenty
  log_prior <-
    dnorm(sigma_treat,   mean = 2, sd = 1, log = TRUE) +
    dnorm(sigma_general, mean = 1, sd = 1, log = TRUE) +
    dnorm(rho_treat,     mean = 3, sd = 1, log = TRUE) +
    dnorm(rho_general,   mean = 0.5, sd = 0.1, log = TRUE)

  list(Score = logML + prior_lambda * log_prior)
}




fit_multiaddgps <- function(sim, hyper_params, seed) {

  sim <- normalize_simulation_design(sim)

  Gamma_general <- function(X) SE(X, hyper_params[3], hyper_params[4], X_row = 3, jitter = 1e-6)
  Gamma_intervention <- function(X) intervention_kernel(
    X, sigma = hyper_params[1], rho = hyper_params[2], X_row = "time", jitter = 1e-6
  )

  linear_component_theta <- matrix(0,sim$D-1,2) 
  linear_component_gamma <- diag(2) + 1e-6 # a matrix with (row and col) dim of the length of linear components

  Theta_kernel <- function(X) matrix(0, sim$D - 1, ncol(X))
  Theta <- list(linear_component_theta, Theta_kernel, Theta_kernel)
  Gamma <- list(linear_component_gamma, Gamma_general, Gamma_intervention)

  start <- Sys.time()
  fit <- fido::basset(
    sim$Y, sim$X,
    upsilon = sim$D + 100,
    Theta = Theta,
    Gamma = Gamma,
    linear = c(1:2),
    n_samples = 2000,
    verbose = TRUE,
    seed = seed + 100,
    jitter = 1e-4,
    max_iter = 1000,
    init = alr_array(sim$Y+0.65, parts = 1)
  )
  runtime <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  fit <- to_clr(fit)

  extract_lambda <- function(idx, Xmult = NULL) {
    arr <- fit$Lambda[idx][[1]]
    if (!is.null(Xmult)) {
      arrr <- arr
      arr <- array(0, dim = c(dim(arr)[1], sim$N, dim(arr)[3]))
      for (i in seq_len(dim(arr)[3])) arr[, , i] <- arrr[, , i] %*% Xmult
    }
    standardize_function(arr, TRUE)
  }

  list(
    estimates = list(
      Cumulative   = tidy_posterior(predict(fit, sim$X, "Lambda",seed = 893), sim),
      Batch        = tidy_posterior(extract_lambda(1, sim$X[1:2, ]), sim),
      General      = tidy_posterior(extract_lambda(2), sim),
      Intervention = tidy_posterior(extract_lambda(3), sim)
    ),
    runtime = runtime
  )
}



fit_multiadd_pipeline <- function(sim, seed) {
  iteration <- 0L
  optimization <- bayesOpt(
    FUN = function(sigma_treat, rho_treat, sigma_general, rho_general) {
      iteration <<- iteration + 1L
      MLL_multiaddgps(
        sim = sim,
        sigma_treat = sigma_treat,
        rho_treat = rho_treat,
        sigma_general = sigma_general,
        rho_general = rho_general,
        seed = seed,
        use_priors = TRUE,
        prior_lambda = 1,
        fast_samplesize = 300,
        fast_max_iter = 200
      )
    },
    bounds = list(
      sigma_treat = c(1, 10),
      rho_treat = c(0.5, 10),
      sigma_general = c(1, 10),
      rho_general = c(0.5, 10)
    ),
    initPoints = 10,
    iters.n = 20,
    verbose = 1
  )
  best <- getBestPars(optimization)
  parameters <- c(
    best$sigma_treat, best$rho_treat,
    best$sigma_general, best$rho_general
  )
  fit <- fit_multiaddgps(sim, parameters, seed)
  list(
    estimates = fit$estimates,
    runtime = fit$runtime,
    diagnostics = list(bo_parameters = parameters)
  )
}
