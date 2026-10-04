# LGPR fitting, prediction summaries, and diagnostics for Simulation 1.

# MCMC diagnostics.
extract_lgpr_diagnostics <- function(fit) {
  stanfit <- get_stanfit(fit)
  summary <- rstan::summary(stanfit)$summary
  parameters <- summary[rownames(summary) != "lp__", , drop = FALSE]

  sampler <- rstan::get_sampler_params(stanfit, inc_warmup = FALSE)
  divergences <- sum(vapply(
    sampler,
    function(x) sum(x[, "divergent__"]),
    numeric(1)
  ))

  list(
    total_postwarmup_draws =
      sum(stanfit@sim$n_save) - sum(stanfit@sim$warmup2),
    min_n_eff = min(parameters[, "n_eff"], na.rm = TRUE),
    median_n_eff = median(parameters[, "n_eff"], na.rm = TRUE),
    mean_n_eff = mean(parameters[, "n_eff"], na.rm = TRUE),
    max_n_eff = max(parameters[, "n_eff"], na.rm = TRUE),
    max_rhat = max(parameters[, "Rhat"], na.rm = TRUE),
    divergences = divergences
  )
}

# Flag problematic LGPR fits.
lgpr_is_flagged <- function(x) {
  !is.finite(x$min_n_eff) || x$median_n_eff < ess_threshold ||
    !is.finite(x$max_rhat) || x$max_rhat > rhat_threshold || x$divergences > 0
}


# Prediction scales used for comparison.
lgpr_summary_spec <- function(data, model, dataset = NULL) {
  spec <- list(
    taxon = 2,
    cumulative_scale = "CLR",
    cumulative_truth = "latent_taxon_2_CLR",
    component_scale = if (model == "lgpr_gaussian")
      "centered_CLR" else "centered_log_count_effect",
    component_truth = "centered_latent_taxon_2_CLR"
  )
  if (model == "lgpr_gaussian") {
    spec$response_mean <- mean(data$clr_Y)
    spec$response_sd <- sd(data$clr_Y)

  } else if (model == "lgpr_nb") {
    other_taxa <- setdiff(seq_len(dataset$D), 2)

    spec$log_count_offset <- log(mean(data$Y)) # baseline

    spec$clr_coefficient <- (dataset$D - 1) / dataset$D
    spec$other_taxa_log_count_mean <-
      colSums(log(dataset$Y[other_taxa, , drop = FALSE] + 0.65)) /
      dataset$D
    spec$pseudocount <- 0.65
    spec$interval_other_taxa <- "observed_and_fixed"

  } else {
    stop("Unsupported LGPR model: ", model)
  }
  spec
}


# Put LGPR posterior draws on the comparison scale.
lgpr_comparison_arrays <- function(arrays, data, model, dataset = NULL) {
  spec <- lgpr_summary_spec(data, model, dataset)

  if (model == "lgpr_gaussian") {
    arrays$Cumulative <-
      spec$response_mean + spec$response_sd * arrays$Cumulative
    for (component in c("Batch", "General", "Intervention")) {
      arrays[[component]] <- spec$response_sd * arrays[[component]]
    }

  } else {
    log_count <- spec$log_count_offset + arrays$Cumulative # add baseline back the the latent effect f
    # Convert the predicted log count for taxon 2 to CLR using the observed counts of the other taxa.
    # CLR_2 = ((D - 1) / D) * log(count_2) - mean log-count contribution from the other taxa.
    arrays$Cumulative <- sweep(spec$clr_coefficient * log_count, 2, spec$other_taxa_log_count_mean,"-")
    
    # Convert centered log-count components to centered CLR effects.
    for (component in c("Batch", "General", "Intervention")) {
      arrays[[component]] <-
        spec$clr_coefficient * arrays[[component]]
    } 
  }
  arrays
}


# Posterior prediction summaries.
extract_lgpr_estimates <- function(
  fit, data, model, draws, seed, dataset = NULL
) {
  stanfit <- get_stanfit(fit)
  available <- sum(stanfit@sim$n_save) - sum(stanfit@sim$warmup2)
  draws <- min(draws, available)

  set.seed(seed + 200)
  chosen <- sample.int(available, draws)
  blocks <- split(chosen, ceiling(seq_along(chosen) / 100))

  get_predictions <- function(indices) {
    invisible(capture.output(
      prediction <- pred(
        fit,
        data,
        draws = indices,
        verbose = FALSE
      )
    ))

    if (model == "lgpr_gaussian") {
      list(
        Cumulative = prediction@f_mean,
        Batch = prediction@f_comp_mean[[2]],
        General = prediction@f_comp_mean[[1]],
        Intervention = prediction@f_comp_mean[[3]]
      )
    } else {
      list(
        Cumulative = prediction@f,
        Batch = prediction@f_comp[[2]],
        General = prediction@f_comp[[1]],
        Intervention = prediction@f_comp[[3]]
      )
    }
  }

  prediction_blocks <- lapply(blocks, get_predictions)

  components <- c("Cumulative", "Batch", "General", "Intervention")
  arrays <- setNames(
    lapply(
      components,
      function(component) {
        do.call(
          rbind,
          lapply(prediction_blocks, `[[`, component)
        )
      }
    ),
    components
  )
  arrays <- lgpr_comparison_arrays(arrays, data, model, dataset)

  tidy_draws <- function(x, center = TRUE) {
    if (center) {
      x <- sweep(x, 1, rowMeans(x), "-")
    }
    q <- apply(x, 2, quantile, probs = c(.025, .25, .5, .75, .975), type = 7)
    data.frame(
      Day = data$time,
      Batch = as.character(data$Batch),
      p2.5 = q[1, ],
      p25 = q[2, ],
      p50 = q[3, ],
      mean = colMeans(x),
      p75 = q[4, ],
      p97.5 = q[5, ]
    )
  }

  list(
    Cumulative = tidy_draws(arrays$Cumulative, center = FALSE),
    Batch = tidy_draws(arrays$Batch),
    General = tidy_draws(arrays$General),
    Intervention = tidy_draws(arrays$Intervention)
  )
}


# Fit LGPR.
fit_lgpr_pipeline <- function(sim, model, iter, chains, cores, seed) {
  if (model == "lgpr_gaussian") {
    data <- make_gam_lgpr_Gaussian_data(sim, taxon_idx = 2)
    formula <- clr_Y ~ gp(Date) + zs(Batch) + gp_ns(Date)
    likelihood <- "gaussian"

  } else if (model == "lgpr_nb") {
    data <- make_gam_lgpr_nb_data(sim, taxon_idx = 2)
    formula <- Y ~ gp(Date) + zs(Batch) + gp_ns(Date)
    likelihood <- "nb"

  } else {
    stop("Unsupported LGPR model: ", model)
  }

  adapt_delta <- as.numeric(
    Sys.getenv("LGPR_ADAPT_DELTA", "0.999")
  )
  max_treedepth <- as.integer(
    Sys.getenv("LGPR_MAX_TREEDEPTH", "11")
  )

  start <- Sys.time()

  fit <- lgp(
    formula,
    data = data,
    likelihood = likelihood,
    prior = list(
      alpha = normal(2, 1),
      ell = normal(0.8, 1),
      wrp = log_normal(1, 1)
    ),
    iter = iter,
    chains = chains,
    cores = cores,
    seed = seed,
    refresh = max(1, floor(iter / 2)),
    verbose = TRUE,
    control = list(
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    )
  )

  runtime <- as.numeric(
    difftime(Sys.time(), start, units = "secs")
  )

  list(
    estimates = extract_lgpr_estimates(
      fit,
      data,
      model,
      draws = 2000,
      seed = seed,
      dataset = sim
    ),
    summary_spec = lgpr_summary_spec(data, model, sim),
    runtime = runtime,
    diagnostics = extract_lgpr_diagnostics(fit),
    sampler_control = list(
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    ),
    fit = fit
  )
}