# Human longitudinal microbiome analysis from David et al. (2014).


library(phyloseq)
library(ape)
library(philr)
library(compositions)
library(ggplot2)
library(ggtree)
library(gridExtra)
library(fido)
library(ParBayesianOptimization)



analysis_dir <-"."
figure_dir <- file.path(analysis_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(analysis_dir, "Human_gut_analysis_helpers.R"), chdir = TRUE)


# Data preparation ---------------

# Aggregate taxa to the family level
physeq <- readRDS(file.path(analysis_dir, "data", "Data1_David", "data_filtered.rds"))
physeq <- tax_glom(physeq, taxrank = "Family")

# Metadata
metadata <- data.frame(sample_data(physeq))
metadata$COLLECTION_DAY <- as.numeric(as.character(metadata$COLLECTION_DAY))
metadata$Day <- as.numeric(scale(metadata$COLLECTION_DAY))

# Design matrix
X <- t(model.matrix(~ Day + COLLECTION_DAY, data = metadata))

# Count matrix
Y <- as(otu_table(physeq), "matrix")
rownames(Y) <- data.frame(tax_table(physeq))$Family

# Dimensions
stopifnot(identical(colnames(Y), colnames(X)))

D <- nrow(Y)
N <- ncol(Y)

cat("Analyzing", D, "families and", N, "samples.\n")
rm(physeq)
gc()


# Multinomial logistic-normal model --------------------------------------

upsilon <- D + 100
G <- cbind(diag(D - 1), -1)
Xi <- (upsilon - D) * base::'%*%'(G, t(G))

fit_model <- function(p, clr = TRUE) {
  zero_mean <- function(X_new) matrix(0, D - 1, ncol(X_new))
  travel_kernel <- function(X_new)
    intervention_kernel(
      X_new, sigma = p$sigma_travel, rho = p$rho_travel, X_row = 3,
      reference_days = metadata$COLLECTION_DAY,
      mode = 92, scale = 20, skew = 2)
  general_kernel <- function(X_new)
    SE(X_new, sigma = p$sigma_general, rho = p$rho_general, X_row = 2)

  fit <- fido::basset(
    Y, X, upsilon,
    Theta = list(matrix(0, D - 1, 1), zero_mean, zero_mean),
    Gamma = list(diag(1) + 1e-6, travel_kernel, general_kernel),
    Xi = Xi, linear = c(1), n_samples = 2000,
    verbose = TRUE, seed = 893, jitter = 1e-4, max_iter = 1000,
    init = fido::alr_array(Y + 0.65, parts = 1)
  )
  if (clr) fido::to_clr(fit) else fit
}

# Hyperparameters ---------------------------------------------------------
# Set to TRUE to re-run Bayesian hyperparameter optimization.

use_optimization <- FALSE

if (use_optimization) {
  score <- function(sigma_travel, rho_travel, sigma_general, rho_general) {
    candidate <- list(
      sigma_travel = sigma_travel,
      rho_travel = rho_travel,
      sigma_general = sigma_general,
      rho_general = rho_general
    )
    constrained <- list(
      sigma_travel = min(sigma_travel, sigma_general),
      rho_travel = max(rho_travel, rho_general),
      sigma_general = max(sigma_travel, sigma_general),
      rho_general = min(rho_travel, rho_general)
    )

    candidate_fit <- fit_model(constrained, clr = FALSE)
    log_prior <- with(candidate,
      dnorm(sigma_travel, 5, 1, log = TRUE) +
      dnorm(sigma_general, 5, 1, log = TRUE) +
      dnorm(rho_travel, 5, 1, log = TRUE) +
      dgamma(1 / rho_general, shape = 10, rate = 8.8, log = TRUE) -
      log(rho_general^2)
    )

    list(Score = candidate_fit$logMarginalLikelihood +
           100 * log_prior)
  }

  set.seed(893)

  optimization <- ParBayesianOptimization::bayesOpt(
    FUN = score,
    bounds = list(
      sigma_travel = c(2, 10),
      rho_travel = c(2, 6),
      sigma_general = c(1, 5),
      rho_general = c(0.5, 2)
    ),
    initPoints = 10,
    iters.n = 20,
    verbose = 1
  )

  best <- ParBayesianOptimization::getBestPars(optimization)
  parameters <- list(
    sigma_travel = min(best$sigma_travel, best$sigma_general),
    rho_travel = max(best$rho_travel, best$rho_general),
    sigma_general = max(best$sigma_travel, best$sigma_general),
    rho_general = min(best$rho_travel, best$rho_general)
  )
} else {
  parameters <- list(
    sigma_travel = 4.02057305496965,
    rho_travel = 4.54129491899598,
    sigma_general = 4.70576775364282,
    rho_general = 0.883663832541547
  )
} # optimized hyperparameters

print(parameters)

# Fit and rank balances ---------------------------------------------------
set.seed(893)
fit <- fit_model(parameters)
analysis <- analyze_balances(fit, metadata, reference_seed = 20261048)

print(analysis$ranking[1:2, c("balance", "percent_mean", "percent_lower","percent_upper", "positive", "negative")],
  row.names = FALSE)

# Figures -----------------------------------------------------------------

tree_plot <- plot_balance_tree(analysis)
ggsave(file.path(figure_dir, "balance_tree.png"), tree_plot,
       width = 12, height = 18, dpi = 300)

top2_plots <- plot_top_balances(analysis, analysis$contrast_matrix,
                  file.path(figure_dir, "top.png"))

# Sensitivity analysis ---------------------------------------------------

run_sensitivity <- TRUE

if (run_sensitivity) {
  # Fix top two balances from the baseline fit
  selected <- analysis$ranking$balance[1:2]
  basis <- analysis$contrast_matrix[, selected, drop = FALSE]

  # Baseline + ±10% perturbation of each hyperparameter
  scenarios <- list(baseline = parameters)
  for (name in names(parameters)) {
    for (mult in c(0.9, 1.1)) {
      p <- parameters
      p[[name]] <- p[[name]] * mult

      suffix <- if (mult == 0.9) "low" else "high"
      scenarios[[paste(name, suffix, sep = "_")]] <- p
    }
  }

  # Fit and summarize each scenario
  summary <- do.call(rbind, lapply(names(scenarios), function(scenario) {
    scenario_fit <- if (scenario == "baseline") {
      fit
    } else {
      fit_model(scenarios[[scenario]])
    }
    summarize_sensitivity(
      scenario,
      scenarios[[scenario]],
      fixed_balance_proportions(scenario_fit, basis)
    )
  }))

  # Plot
  labels <- sensitivity_scenario_labels()
  summary$scenario_factor <- factor(
    summary$scenario,
    levels = rev(names(labels))
  )
  summary$is_baseline <- summary$scenario == "baseline"
  for (balance in selected) {
    plot_sensitivity_balance(
      summary, selected, balance, labels, figure_dir
    )
  }
}

# Priors used below: Normal(5, 1) for both amplitudes and travel length scale;
# InverseGamma(10, 8.8) for the general length scale.
dir.create(figure_dir, showWarnings = FALSE)
png(file.path(figure_dir, "Web_Figure9.png"), width = 2100, height = 1200, res = 300)
par(mfrow = c(1, 2))
prior_x <- seq(0, 10, length.out = 1000)
plot(prior_x, dnorm(prior_x, mean = 5, sd = 1),
     type = "l", col = "darkred", lwd = 1.5, xlim = c(0, 10), las = 1,
     xlab = "Normal(5, 1)", ylab = "Prior Density")
prior_x <- seq(0.001, 4, length.out = 1000)
plot(prior_x, dgamma(1 / prior_x, shape = 10, rate = 8.8) / prior_x^2,
     type = "l", col = "darkred", lwd = 1.5, xlim = c(0, 4), las = 1,
     xlab = "InverseGamma(10, 8.8)", ylab = "")
dev.off()
