# Functions used by simulation2_fit.R and simulation2.R.
# This file only defines functions; sourcing it does not fit a model.

# Kernels ---------------------------------------------------------------

SE <- function(X, sigma, rho, X_row = 3, jitter = 1e-10) {
  distance <- as.matrix(dist(t(X[X_row, , drop = FALSE])))
  sigma^2 * exp(-distance^2 / (2 * rho^2)) + jitter * diag(ncol(distance))
}

Period <- function(X, sigma, rho, p = 30, X_row = 3, jitter = 1e-10) {
  distance <- as.matrix(dist(t(X[X_row, , drop = FALSE])))
  sigma^2 * exp(-2 * sin(pi * distance / p)^2 / rho^2) +
    jitter * diag(ncol(distance))
}

# Model matrices and marginal likelihoods -------------------------------

order_parameters <- function(p) {
  list(sigma1 = max(p$sigma1, p$sigma2), l1 = min(p$l1, p$l2),
       sigma2 = min(p$sigma1, p$sigma2), l2 = max(p$l1, p$l2))
}

model_matrices <- function(X, parameters, beta_0, beta_1, Xi, upsilon) {
  p <- order_parameters(parameters)
  q <- nrow(Xi)
  Xb <- X[1:2, , drop = FALSE]
  Bmean <- cbind(rep(beta_0, q), rep(beta_1, q))
  # Keep the original addition to every entry, not only the diagonal.
  Gb <- diag(2) * 0.1 + 1e-6
  Kp <- Period(X, sigma = p$sigma1, rho = p$l1, p = 30)
  Kt <- SE(X, sigma = p$sigma2, rho = p$l2)
  G <- t(Xb) %*% Gb %*% Xb + Kp + Kt
  list(parameters = p, Xb = Xb, Bmean = Bmean, Gb = Gb, Kp = Kp,
       Kt = Kt, Gamma = G, Theta = Bmean %*% Xb,
       Xi = Xi, upsilon = upsilon, A_chol = chol(diag(ncol(X)) + G))
}

log_multigamma <- function(a, q) {
  q * (q - 1) / 4 * log(pi) + sum(lgamma(a + (1 - seq_len(q)) / 2))
}

# Integrate out F and Sigma under the shared inverse-Wishart prior.
matrix_t_log_density <- function(H, model) {
  q <- nrow(H)
  n <- ncol(H)
  residual <- H - model$Theta
  whitened <- forwardsolve(t(model$A_chol), t(residual))
  post_scale <- model$Xi + crossprod(whitened)
  nu <- model$upsilon
  log_multigamma((nu + n) / 2, q) - log_multigamma(nu / 2, q) -
    q * n / 2 * log(pi) - q * sum(log(diag(model$A_chol))) +
    nu * sum(log(diag(chol(model$Xi)))) -
    (nu + n) * sum(log(diag(chol(post_scale))))
}

hyperparameter_log_penalty <- function(parameters, weight = 100) {
  p <- order_parameters(parameters)
  weight * (dnorm(p$sigma1, 4, 1, log = TRUE) +
              dnorm(p$sigma2, 1, 1, log = TRUE) +
              dnorm(p$l1, 25, 1, log = TRUE) +
              dnorm(p$l2, 30, 1, log = TRUE))
}

multinomial_log_density <- function(Y, H) {
  # Include the ALR reference category with logit zero; use log-sum-exp.
  peak <- pmax(0, apply(H, 2, max))
  log_denom <- peak + log(exp(-peak) + colSums(exp(sweep(H, 2, peak))))
  sum(lgamma(colSums(Y) + 1)) - sum(lgamma(Y + 1)) +
    sum(Y[-nrow(Y), , drop = FALSE] * H) - sum(colSums(Y) * log_denom)
}

collapse_multiadd <- function(Y, model, n_samples = 0) {
  # Convert the inverse-Wishart df to fido's matrix-t df (fido 1.1.1).
  native_df <- model$upsilon - (nrow(Y) - 1) + 1
  native <- getFromNamespace("optimPibbleCollapsed", "fido")(
    Y = Y, upsilon = native_df, ThetaX = model$Theta,
    KInv = chol2inv(chol(model$Xi)), AInv = chol2inv(model$A_chol),
    init = fido::alr_array(Y + 0.65, parts = 1),
    # This fido release skips the Hessian determinant at n_samples = 0.
    # One disposable draw obtains the same Laplace score used in the final fit.
    n_samples = as.integer(max(1, n_samples)), calcGradHess = TRUE,
    max_iter = 1000, jitter = 1e-4,
    ncores = 1, seed = 893, verbose = FALSE
  )
  stopifnot(all(is.finite(native$Pars)), length(native$logInvNegHessDet) == 1,
            is.finite(native$logInvNegHessDet))
  # Use the normalized marginal density for the SAME IW prior as NAddGP.
  score <- multinomial_log_density(Y, native$Pars) +
    matrix_t_log_density(native$Pars, model) +
    (nrow(Y) - 1) * ncol(Y) / 2 * log(2 * pi) +
    native$logInvNegHessDet / 2
  if (n_samples > 0) {
    stopifnot(identical(dim(native$Samples), as.integer(c(nrow(Y) - 1, ncol(Y), n_samples))))
  }
  list(log_marginal = score, H = native$Samples, map = native$Pars,
       max_abs_gradient = max(abs(native$Gradient)))
}

# Posterior sampling ----------------------------------------------------

sample_F <- function(H, model, seed = 894) {
  # Sigma | H ~ IW(Xi + R A^{-1} R^T, upsilon + N).
  # F | H, Sigma ~ MN(H - R A^{-1}, Sigma, I - A^{-1}).
  set.seed(seed)
  q <- dim(H)[1]
  n <- dim(H)[2]
  S <- dim(H)[3]
  Ainv <- chol2inv(model$A_chol)
  C <- diag(n) - Ainv
  root_C <- chol((C + t(C)) / 2)
  F <- array(NA_real_, dim(H))
  Sigma <- array(NA_real_, c(q, q, S))
  for (s in seq_len(S)) {
    residual <- H[, , s] - model$Theta
    residual_Ainv <- residual %*% Ainv
    post_scale <- model$Xi + residual_Ainv %*% t(residual)
    post_scale <- (post_scale + t(post_scale)) / 2
    precision <- rWishart(1, model$upsilon + n,
                          chol2inv(chol(post_scale)))[, , 1]
    Sigma[, , s] <- chol2inv(chol(precision))
    F[, , s] <- H[, , s] - residual_Ainv +
      t(chol(Sigma[, , s])) %*% matrix(rnorm(q * n), q) %*% root_C
  }
  list(Lambda = F, Sigma = Sigma)
}

sample_components <- function(F, Sigma, model, seed = 895) {
  # Draw independent prior components, then condition their sum on F.
  # Right multiplication by chol(K) gives column covariance K.
  set.seed(seed)
  Gchol <- chol(model$Gamma)
  solve_G <- function(rhs) backsolve(Gchol, forwardsolve(t(Gchol), rhs))
  gain_B <- solve_G(t(model$Xb) %*% model$Gb)
  gain_P <- solve_G(model$Kp)
  root_B <- chol(model$Gb)
  root_P <- chol(model$Kp)
  root_T <- chol(model$Kt)
  q <- dim(F)[1]
  n <- dim(F)[2]
  S <- dim(F)[3]
  B <- array(NA_real_, c(q, ncol(model$Bmean), S))
  batch <- periodic <- trend <- array(NA_real_, dim(F))
  for (s in seq_len(S)) {
    row_root <- t(chol(Sigma[, , s]))
    prior_B <- model$Bmean + row_root %*% matrix(rnorm(q * 2), q) %*% root_B
    prior_P <- row_root %*% matrix(rnorm(q * n), q) %*% root_P
    prior_T <- row_root %*% matrix(rnorm(q * n), q) %*% root_T
    residual <- F[, , s] - prior_B %*% model$Xb - prior_P - prior_T
    B[, , s] <- prior_B + residual %*% gain_B
    batch[, , s] <- B[, , s] %*% model$Xb
    periodic[, , s] <- prior_P + residual %*% gain_P
    trend[, , s] <- F[, , s] - batch[, , s] - periodic[, , s]
  }
  list(B = B, components = list(batch, periodic, trend, F),
       sum_error = max(abs(batch + periodic + trend - F)))
}

# Posterior summaries and plotting -------------------------------------

center_over_time <- function(x) {
  if (length(dim(x)) == 3) {
    sweep(x, c(1, 3), apply(x, c(1, 3), mean), "-")
  } else {
    sweep(x, 1, rowMeans(x), "-")
  }
}

summarize_coordinate <- function(x, coordinate, X1) {
  # All archived time points are unique. Summarize the selected coordinate
  # directly, avoiding a multi-million-row intermediate table.
  stopifnot(!anyDuplicated(X1$time))
  if (length(dim(x)) == 3) {
    draws <- x[coordinate, , ]
    quantiles <- t(apply(draws, 1, quantile, c(.025, .25, .5, .75, .975)))
    out <- data.frame(mean = rowMeans(draws), p2.5 = quantiles[, 1],
                      p25 = quantiles[, 2], p50 = quantiles[, 3],
                      p75 = quantiles[, 4], p97.5 = quantiles[, 5])
  } else {
    values <- x[coordinate, ]
    out <- data.frame(mean = values, p2.5 = values, p25 = values,
                      p50 = values, p75 = values, p97.5 = values)
  }
  out$Date <- X1$time
  out$batch <- X1$batch
  out[order(out$Date, out$batch), ]
}

summarize_component <- function(draws, coordinate, X1, suffix = "") {
  summaries <- lapply(names(draws), function(source) {
    result <- summarize_coordinate(draws[[source]], coordinate, X1)
    result$source <- paste0(source, suffix)
    result
  })

  do.call(rbind, summaries)
}

component_plot <- function(data, title, subtitle, figure_theme, custom_colors,
                           batch_effect = FALSE) {
  if (batch_effect) {
    plot <- ggplot(data, aes(x = Date, y = mean, group = interaction(batch, source))) +
      geom_ribbon(
        aes(ymin = p2.5, ymax = p97.5, fill = source, group = interaction(batch, source)),
        alpha = 0.2
      ) +
      geom_ribbon(
        aes(ymin = p25, ymax = p75, fill = source, group = interaction(batch, source)),
        alpha = 0.5
      ) +
      geom_line(aes(color = source, linetype = batch), linewidth = 1)
  } else {
    plot <- ggplot(data, aes(x = Date, y = mean, group = source)) +
      geom_ribbon(aes(ymin = p2.5, ymax = p97.5, fill = source, group = source), alpha = 0.2) +
      geom_ribbon(aes(ymin = p25, ymax = p75, fill = source, group = source), alpha = 0.5) +
      geom_line(aes(color = source), linewidth = 0.7)
  }

  plot +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = custom_colors) +
    figure_theme +
    theme(legend.position = "none") +
    labs(y = "F", x = "Hour", title = title, subtitle = subtitle)
}
