# Run from the folder containing these three scripts.
# This script fits the models. Use simulation2.R to plot saved results.

library(fido)
library(ParBayesianOptimization)

project_dir <- "."
data_dir <- file.path(project_dir, "data/Data0_Simulation/data")
result_dir <- data_dir
source(file.path(project_dir, "simulation2_helper.R"))
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# Data and priors -------------------------------------------------------

# Load the archived observations; do not regenerate the simulation.
data_env <- new.env()
load(file.path(data_dir, "simulated_data2.RData"), envir = data_env)
sim <- data_env$simulated_data2
Y <- sim$Y
X <- sim$X
D <- as.integer(sim$D)
N <- as.integer(sim$N)
beta_0 <- sim$beta_0
beta_1 <- sim$beta_1
q <- D - 1
samplesize <- 2000
stopifnot(D == 4, N == 600, all(colSums(Y) == 20))

upsilon <- D + 1000
prior_correlation <- matrix(0.9, q, q)
diag(prior_correlation) <- 1
Xi <- (upsilon - D - 1) * prior_correlation
H_fixed <- fido::alr_array(Y + 0.65, parts = 1)
bounds <- list(sigma1 = c(1, 9), l1 = c(15, 35),
               sigma2 = c(1, 6), l2 = c(20, 40))

# NAddGP: optimize the hyperparameters ----------------------------------

score_NAddGP <- function(sigma1, l1, sigma2, l2) {
  parameters <- list(sigma1 = sigma1, l1 = l1, sigma2 = sigma2, l2 = l2)
  model <- model_matrices(X, parameters, beta_0, beta_1, Xi, upsilon)
  list(Score = matrix_t_log_density(H_fixed, model) +
         hyperparameter_log_penalty(parameters))
}

set.seed(893)
optObj_NAddGP <- bayesOpt(
  FUN = score_NAddGP, bounds = bounds,
  initPoints = 10, iters.n = 40, verbose = 1
)
raw_parameters <- as.list(getBestPars(optObj_NAddGP))
params <- order_parameters(raw_parameters)
nadd_checked_score <- do.call(score_NAddGP, params)$Score
stopifnot(is.finite(nadd_checked_score),
          abs(nadd_checked_score - max(optObj_NAddGP$scoreSummary$Score)) < 1e-5)
print(params)

# MultiAddGP: optimize the hyperparameters -------------------------------

score_MultiAddGP <- function(sigma1, l1, sigma2, l2) {
  parameters <- list(sigma1 = sigma1, l1 = l1, sigma2 = sigma2, l2 = l2)
  model <- model_matrices(X, parameters, beta_0, beta_1, Xi, upsilon)
  collapsed <- collapse_multiadd(Y, model)
  list(Score = collapsed$log_marginal + hyperparameter_log_penalty(parameters))
}

set.seed(893)
optObj_MultiAddGP <- bayesOpt(
  FUN = score_MultiAddGP, bounds = bounds,
  initPoints = 10, iters.n = 40, verbose = 1
)
raw_parameters <- as.list(getBestPars(optObj_MultiAddGP))
hyper <- order_parameters(raw_parameters)
multi_checked_score <- do.call(score_MultiAddGP, hyper)$Score
stopifnot(is.finite(multi_checked_score),
          abs(multi_checked_score - max(optObj_MultiAddGP$scoreSummary$Score)) < 1e-5)
print(hyper)


# NAddGP: sample F, Sigma, and the three components -----------------------

model_NAddGP <- model_matrices(X, params, beta_0, beta_1, Xi, upsilon)
H <- array(rep(H_fixed, samplesize), c(q, N, samplesize))
nadd <- sample_F(H, model_NAddGP)
nadd_parts <- sample_components(nadd$Lambda, nadd$Sigma, model_NAddGP)
sample_posterior <- nadd_parts$components

# MultiAddGP: sample H, then F, Sigma, and the three components ------------

model_MultiAddGP <- model_matrices(X, hyper, beta_0, beta_1, Xi, upsilon)
collapsed <- collapse_multiadd(Y, model_MultiAddGP, n_samples = samplesize)
stopifnot(abs(collapsed$log_marginal + hyperparameter_log_penalty(hyper) -
                multi_checked_score) < 1e-5)
multi <- sample_F(collapsed$H, model_MultiAddGP)
multi_parts <- sample_components(multi$Lambda, multi$Sigma, model_MultiAddGP)

mod <- list(
  Lambda = list(multi_parts$B, multi_parts$components[[2]],
                multi_parts$components[[3]])
)
predicted <- multi_parts$components[[4]]
stopifnot(nadd_parts$sum_error < 1e-9, multi_parts$sum_error < 1e-9)

# Save only the posterior draws needed for Web Figures 3-5 --------------

saveRDS(
  list(mod = mod, sample_posterior = sample_posterior,
       predicted = predicted),
  file.path(result_dir, "Simulation_result2.rds")
)
