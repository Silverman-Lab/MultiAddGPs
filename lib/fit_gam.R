# GAM fitting and component summaries for Simulation 1.

tidy_gam_summary <- function(data, mean, se) {
  data.frame(
    Day = data$time,
    Batch = as.character(data$Batch),
    p2.5 = mean - 1.96 * se,
    p25 = mean - 0.674 * se,
    p50 = mean,
    mean = mean,
    p75 = mean + 0.674 * se,
    p97.5 = mean + 1.96 * se
  )
}

centered_gam_summary <- function(data, design, coefficients, covariance) {
  centered_design <- sweep(design, 2L, colMeans(design), "-")
  variance <- rowSums((centered_design %*% covariance) * centered_design)
  stopifnot(all(is.finite(variance)), min(variance) >= -1e-8)
  tidy_gam_summary(
    data,
    drop(centered_design %*% coefficients),
    sqrt(pmax(variance, 0))
  )
}

summarize_gam_fit <- function(gam_model, data) {
  total <- predict(gam_model, se.fit = TRUE)
  design <- predict(gam_model, type = "lpmatrix")
  coefficients <- coef(gam_model)
  covariance <- gam_model$Vp

  batch_data <- data.frame(
    Date = 0, Date2 = 0, Batch = c("Batch1", "Batch2")
  )
  batch_index <- match(as.character(data$Batch), batch_data$Batch)
  stopifnot(!anyNA(batch_index))
  batch_design <- predict(
    gam_model, newdata = batch_data, type = "lpmatrix"
  )[batch_index, , drop = FALSE]

  smooth_summary <- function(label) {
    index <- which(vapply(
      gam_model$smooth, function(smooth) identical(smooth$label, label),
      logical(1)
    ))
    stopifnot(length(index) == 1L)
    smooth <- gam_model$smooth[[index]]
    columns <- seq.int(smooth$first.para, smooth$last.para)
    centered_gam_summary(
      data, design[, columns, drop = FALSE], coefficients[columns],
      covariance[columns, columns, drop = FALSE]
    )
  }

  list(
    Cumulative = tidy_gam_summary(data, total$fit, total$se.fit),
    Batch = centered_gam_summary(
      data, batch_design, coefficients, covariance
    ),
    General = smooth_summary("s(Date)"),
    Intervention = smooth_summary("s(Date2)")
  )
}

fit_gam <- function(data, time) {
  start <- Sys.time()
  gam_model <- gam(
    clr_Y ~ Batch + s(Date, k = 5) + s(Date2, bs = "ad", k = 12),
    data = data
  )
  runtime <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  list(
    estimates = summarize_gam_fit(gam_model, data),
    runtime = runtime
  )
}
