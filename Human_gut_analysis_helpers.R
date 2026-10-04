
# The squared-exponential kernel
SE <- function(X, sigma = 1, rho = 1, X_row, jitter = 1e-6) {
  distance <- as.matrix(stats::dist(t(X[X_row, , drop = FALSE])))
  sigma^2 * exp(-distance^2 / (2 * rho^2)) + jitter * diag(ncol(X))
}

# Intervention time transformation.
# Normalize using a fixed reference set of observed collection days.
# Keep reference_days unchanged when evaluating subsets or new dates.
warp <- function(x, reference_days, mode = 92, scale = 20, skew = 2) {
  unscaled <- function(t) {
    z <- (t - mode) / scale
    exp(-0.5 * z^2) * (1 / (1 + exp(-skew * z)))
  }
  limits <- range(unscaled(reference_days))
  2 * (unscaled(x) - limits[1]) / diff(limits) - 1
}

# Intervention covariance from raw collection days.
# X has samples in columns; X_row selects the raw collection-day row.
# reference_days is the full observed collection-day vector for normalization.
intervention_kernel <- function(X, sigma = 1, rho = 1,
                                X_row, reference_days, jitter = 1e-6,
                                mode = 92, scale = 20, skew = 2) {
  X[X_row, ] <- warp(X[X_row, ], reference_days = reference_days,
                    mode = mode, scale = scale, skew = skew)
  SE(X, sigma = sigma, rho = rho, X_row = X_row, jitter = jitter)
}

centered_travel <- function(fit, s) {
  # Row centering removes the draw-specific intercept from this component.
  L <- fit$Lambda[[2]][, , s] + fit$Lambda[[1]][, , s]
  L - rowMeans(L)
}

build_ward_trees <- function(fit, taxa) {
  S <- dim(fit$Lambda[[2]])[3]
  trees <- vector("list", S)
  for (s in seq_len(S)) {
    W <- stats::cov(t(centered_travel(fit, s)))
    # Pairwise variance of taxon differences is the Ward-tree dissimilarity.
    variation <- outer(diag(W), diag(W), "+") - 2 * W
    dimnames(variation) <- NULL
    trees[[s]] <- ape::as.phylo(stats::hclust(stats::as.dist(variation), "ward.D"))
    trees[[s]]$tip.label <- taxa[as.integer(trees[[s]]$tip.label)]
  }
  trees
}

make_consensus_balance_basis <- function(trees, taxa, reference_seed) {
  consensus <- ape::consensus(trees, p = 0.5)
  # RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(reference_seed)
  tree <- ape::makeNodeLabel(ape::multi2di(consensus, random = TRUE),
                             method = "number", prefix = "w")
  V <- compositions::gsi.buildilrBase(philr::phylo2sbp(tree))
  colnames(V) <- tree$node.label
  V <- V[taxa, , drop = FALSE]
  list(tree = tree, consensus_internal_nodes = consensus$Nnode,
       contrast_matrix = V)
}

calculate_balance_variation <- function(fit, V) {
  D <- nrow(fit$Y)
  S <- dim(fit$Lambda[[2]])[3]
  variances <- matrix(NA_real_, D - 1L, S, dimnames = list(colnames(V), NULL))
  total <- numeric(S)
  for (s in seq_len(S)) {
    L <- centered_travel(fit, s)
    variances[, s] <- apply(crossprod(V, L), 1, stats::var)
    total[s] <- sum(apply(L, 1, stats::var))
  }
  proportions <- sweep(variances, 2, total, "/")
  list(variances = variances, proportions = proportions)
}

rank_balances <- function(variances, proportions, V, taxa) {
  ord <- order(rowMeans(variances), decreasing = TRUE)
  ci <- t(apply(proportions, 1, quantile, c(0.025, 0.975)))
  data.frame(
    rank = seq_along(ord), balance = colnames(V)[ord],
    mean_temporal_variance = rowMeans(variances)[ord],
    percent_mean = 100 * rowMeans(proportions)[ord],
    percent_lower = 100 * ci[ord, 1], percent_upper = 100 * ci[ord, 2],
    positive = vapply(ord, function(k) paste(taxa[V[, k] > 0], collapse = "; "), ""),
    negative = vapply(ord, function(k) paste(taxa[V[, k] < 0], collapse = "; "), ""),
    row.names = NULL
  )
}

calculate_balance_trajectories <- function(fit, metadata, V, selected) {
  S <- dim(fit$Lambda[[2]])[3]
  draws <- array(NA_real_, c(2L, ncol(fit$Y), S))
  for (s in seq_len(S)) {
    draws[, , s] <- crossprod(V[, selected, drop = FALSE], centered_travel(fit, s))
  }
  trajectories <- do.call(rbind, lapply(seq_along(selected), function(k) {
    b <- draws[k, , ]
    q <- t(apply(b, 1, quantile, c(0.025, 0.25, 0.75, 0.975)))
    data.frame(balance = selected[k], day = metadata$COLLECTION_DAY,
               mean = rowMeans(b), lower = q[, 1], q25 = q[, 2],
               q75 = q[, 3], upper = q[, 4])
  }))
  trajectories$balance <- factor(trajectories$balance, levels = selected)
  trajectories <- trajectories[order(trajectories$balance, trajectories$day), ]
  trajectories
}

analyze_balances <- function(fit, metadata, reference_seed) {
  taxa <- rownames(fit$Y)

  # 1. Cluster families separately in each posterior draw.
  trees <- build_ward_trees(fit, taxa)

  # 2. Resolve the consensus tree and turn it into balance directions.
  basis <- make_consensus_balance_basis(trees, taxa, reference_seed)
  V <- basis$contrast_matrix

  # 3. Measure each balance's share of temporal variation in every draw.
  variation <- calculate_balance_variation(fit, V)
  proportions <- variation$proportions

  # 4. Rank balances by mean temporal variance and select the top two.
  ranking <- rank_balances(variation$variances, proportions, V, taxa)
  selected <- ranking$balance[1:2]

  # 5. Summarize the selected balance trajectories.
  trajectories <- calculate_balance_trajectories(fit, metadata, V, selected)
  list(tree = basis$tree,
       consensus_internal_nodes = basis$consensus_internal_nodes,
       contrast_matrix = V, ranking = ranking, proportion_draws = proportions,
       trajectories = trajectories)
}

plot_balance_tree <- function(analysis) {
  p <- ggtree(analysis$tree, branch.length = "none", ladderize = FALSE,
              linewidth = 0.7)
  depth <- max(p$data$x)
  p + geom_tiplab(geom = "label", size = 4, fill = "white") +
    geom_label(data = p$data[!p$data$isTip, ], aes(x = x, y = y, label = label),
               size = 3, fill = "white") +
    xlim(-1, depth * 1.9) + scale_y_reverse() + theme_tree() +
    theme(plot.margin = margin(10, 20, 10, 10))
}

# This display schematic shows the two sides of a balance. The resolved
# consensus tree, with its full topology, is shown in balance_tree.png.
balance_display_tree <- function(V, balance, max_named_families = 5L) {
  taxa <- rownames(V)
  negative <- taxa[V[, balance] < 0]
  positive <- taxa[V[, balance] > 0]
  display <- function(families) {
    if (length(families) <= max_named_families) families
    else sprintf("Other %d families", length(families))
  }
  # Tree tips are drawn in reverse order on the displayed y-axis.
  negative_labels <- rev(display(negative))
  positive_labels <- rev(display(positive))
  labels <- c(negative_labels, positive_labels)
  ids <- paste0("t", seq_along(labels))
  clade <- function(tips) {
    if (length(tips) == 1L) tips else paste0("(", paste(tips, collapse = ","), ")")
  }
  newick <- sprintf("(%s,%s);",
                    clade(ids[seq_along(negative_labels)]),
                    clade(ids[length(negative_labels) + seq_along(positive_labels)]))
  tree <- ape::read.tree(text = newick)
  tree$tip.label <- labels
  list(tree = tree,
       negative_y = mean(seq_along(negative_labels)),
       positive_y = length(negative_labels) + mean(seq_along(positive_labels)))
}

balance_tree_panel <- function(diagram) {
  tree <- diagram$tree
  n_tips <- length(tree$tip.label)
  x_limits <- if (n_tips == 2L) c(-2, 6) else c(-1.5, 7)
  y_limits <- if (n_tips == 2L) c(-0.15, 3.15) else c(0, n_tips + 0.5)
  ggtree::ggtree(tree, linewidth = 1.5, ladderize = FALSE) +
    ggtree::geom_tiplab(geom = "label", align = FALSE, fill = "white",
                        size = 4, offset = 0) +
    xlim(x_limits[1], x_limits[2]) + ggtree::theme_tree2() +
    theme(plot.margin = margin(5, 0, 20, 0),
          axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), axis.line.x = element_blank()) +
    annotate("segment", x = -0.5, xend = -0.5,
             y = diagram$negative_y, yend = diagram$positive_y,
             arrow = grid::arrow(length = grid::unit(0.25, "cm"), ends = "both"),
             linewidth = 0.8, color = "gray40") +
    annotate("text", x = -1, y = diagram$positive_y, label = "+", size = 6,
             fontface = "bold", color = "gray40") +
    annotate("text", x = -1, y = diagram$negative_y, label = "-", size = 6,
             fontface = "bold", color = "gray40") +
    annotate("text", x = -1,
             y = (diagram$negative_y + diagram$positive_y) / 2,
             label = "Balance", angle = 90, size = 3.3,
             color = "gray40") +
    coord_cartesian(ylim = y_limits, clip = "off") + scale_y_reverse()
}

balance_trajectory_panel <- function(analysis, ranking) {
  trajectories <- analysis$trajectories
  labels <- setNames(sprintf("%s (variance rank %d; %.1f%% explained)",
                             ranking$balance, ranking$rank,
                             ranking$percent_mean), ranking$balance)
  trajectories$panel <- factor(labels[as.character(trajectories$balance)],
                               levels = labels[ranking$balance])
  trajectory_plot <- ggplot(trajectories, aes(day, mean)) +
    facet_wrap(~panel, ncol = 1) +
    geom_rect(inherit.aes = FALSE, xmin = 71, xmax = 122,
              ymin = -Inf, ymax = Inf, fill = "#f3efeff1", alpha = 0.1) +
    annotate("segment", x = 71, xend = 122, y = 0.2, yend = 0.2,
             arrow = grid::arrow(ends = "both", length = grid::unit(0.2, "cm")),
             linewidth = 0.4, color = "black") +
    annotate("text", x = 96.5, y = 0.2, label = "Travel", vjust = -0.6,
             size = 3, fontface = "italic") +
    coord_cartesian(clip = "off") +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#2a9cf8", alpha = 0.3) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#2a9cf8", alpha = 0.8) +
    geom_line(color = "#eb176f", linewidth = 0.8) +
    scale_x_continuous(breaks = c(0, 100, 200, 300),
                       labels = function(x) paste0("Day ", x)) +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_line(color = "gray90"),
          panel.grid.minor = element_line(color = "gray95"),
          strip.background = element_rect(fill = "#d4d4d4", color = "black"),
          strip.text = element_text(color = "black", size = 13),
          panel.spacing = grid::unit(1, "lines"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          axis.title.y = element_text(size = 12, margin = margin(t = 15)),
          axis.title.x = element_text(size = 15, margin = margin(t = 15)),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
          axis.text.y = element_text(angle = 0, hjust = 1, size = 10),
          plot.title.position = "plot",
          plot.title = element_text(size = 17, face = "bold")) +
    labs(x = NULL, y = " Balance Value")
  trajectory_plot
}

plot_top_balances <- function(analysis, V, outfile) {
  selected <- analysis$ranking$balance[1:2]
  diagrams <- lapply(selected, function(balance) balance_display_tree(V, balance))
  panels <- lapply(diagrams, balance_tree_panel)
  trajectory_plot <- balance_trajectory_panel(analysis, analysis$ranking[1:2, ])

  plot <- gridExtra::arrangeGrob(panels[[1]], panels[[2]], trajectory_plot,
                                 layout_matrix = rbind(c(1, 3), c(2, 3)),
                                 widths = c(0.6, 0.8))
  grDevices::png(outfile, width = 8, height = 6, units = "in", res = 300,
                 type = "cairo")
  tryCatch({
    grid::grid.newpage()
    grid::grid.draw(plot)
    grid::grid.text("A", x = 0.025, y = 0.985, just = c("left", "top"),
                    gp = grid::gpar(fontsize = 17, fontface = "bold"))
    grid::grid.text("B", x = 0.445, y = 0.985, just = c("left", "top"),
                    gp = grid::gpar(fontsize = 17, fontface = "bold"))
  }, finally = grDevices::dev.off())
  message("Saved ", outfile)
}

fixed_balance_proportions <- function(fit, basis) {
  draws <- matrix(NA_real_, dim(fit$Lambda[[2]])[3], 3L,
                  dimnames = list(NULL, c(colnames(basis), "combined")))
  for (s in seq_len(nrow(draws))) {
    L <- centered_travel(fit, s)
    total <- sum(apply(L, 1L, stats::var))
    projected <- crossprod(basis, L)
    draws[s, 1:2] <- apply(projected, 1L, stats::var) / total
    draws[s, 3L] <- sum(draws[s, 1:2])
  }
  draws
}

summarize_sensitivity <- function(scenario, parameters, draws) {
  data.frame(
    scenario = scenario,
    balance = colnames(draws),
    mean_percent = 100 * colMeans(draws),
    lower_95 = 100 * apply(draws, 2L, stats::quantile, 0.025),
    upper_95 = 100 * apply(draws, 2L, stats::quantile, 0.975),
    sigma_travel = parameters$sigma_travel,
    rho_travel = parameters$rho_travel,
    sigma_general = parameters$sigma_general,
    rho_general = parameters$rho_general,
    row.names = NULL
  )
}

sensitivity_scenario_labels <- function() {
  c(baseline = "Baseline",
    sigma_travel_low = "σ travel: -10%",
    sigma_travel_high = "σ travel: +10%",
    rho_travel_low = "ρ travel: -10%",
    rho_travel_high = "ρ travel: +10%",
    sigma_general_low = "σ general: -10%",
    sigma_general_high = "σ general: +10%",
    rho_general_low = "ρ general: -10%",
    rho_general_high = "ρ general: +10%")
}

plot_sensitivity_balance <- function(summary, selected, balance, labels, final_dir) {
  plot_data <- summary[summary$balance == balance, ]
  baseline_mean <- plot_data$mean_percent[plot_data$scenario == "baseline"]
  plot_title <- paste("Sensitivity analysis for fixed balance", balance)
  plot <- ggplot(plot_data, aes(y = scenario_factor)) +
    geom_vline(xintercept = baseline_mean, linetype = "dashed",
               linewidth = 0.45, color = "gray45") +
    geom_segment(aes(x = lower_95, xend = upper_95,
                     yend = scenario_factor),
                 linewidth = 0.8, color = "gray35") +
    geom_point(aes(x = mean_percent, fill = is_baseline),
               shape = 21, size = 2.8, stroke = 0.65, color = "black") +
    geom_text(aes(x = upper_95,
                  label = sprintf("%.1f%%", mean_percent)),
              hjust = -0.15, size = 3.1, color = "gray20") +
    scale_y_discrete(labels = function(x) unname(labels[x]),
                     expand = expansion(mult = c(0.02, 0.02))) +
    scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white"),
                      guide = "none") +
    labs(title = plot_title,
         subtitle = paste("Points show posterior means; horizontal bars",
                          "show 95% credible intervals"),
         x = "Percentage of total variation explained", y = NULL) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9, color = "gray35"),
          axis.text.y = element_text(size = 9.5, color = "black"),
          axis.text.x = element_text(size = 9.5, color = "black"),
          axis.title.x = element_text(size = 10.5,
                                      margin = margin(t = 7)),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_line(color = "gray90",
                                            linewidth = 0.35),
          plot.margin = margin(5, 30, 5, 5)) +
    coord_cartesian(clip = "off")
  ggsave(file.path(final_dir, paste0("sensitivity_", balance, ".png")),
         plot, width = if (identical(balance, selected[2])) 7 else 8,
         height = 4.5, dpi = 300)
}
