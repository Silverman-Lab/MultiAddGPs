# Truth and posterior formatting for Simulation 1.
#
# These helpers convert simulated truth and posterior arrays into the common
# tidy summaries used by the comparison plots and metric calculations.

extract_truth <- function(sim, coord = "2") {

  tidy_truth <- function(arr) {
    N <- ncol(arr)
    gather_array(arr, mean, coord, sample) %>%
      mutate(
        Day = seq(0,N-1,by = 1)[sample] ,
        Batch = ifelse(sim$X[2, sample] == 1, "Batch1", "Batch2")
      ) %>%
      filter(coord == !!coord) %>%
      group_by(Day, Batch) %>%
      summarise_posterior(mean, na.rm = TRUE) %>%
      ungroup() %>%
      dplyr::select(Day, Batch, mean)
  }

  list(
    Cumulative   = tidy_truth(alr_to_clr(sim$F, sim$D)),
    Batch        = tidy_truth(standardize_function(alr_to_clr(sim$f1X1, sim$D), FALSE)),
    General      = tidy_truth(standardize_function(alr_to_clr(sim$f2, sim$D), FALSE)),
    Intervention = tidy_truth(standardize_function(alr_to_clr(sim$f3, sim$D), FALSE))
  )
}

tidy_posterior <- function(arr, sim, coord_keep = "2") {

  gather_array(arr, val, coord, sample, iter) %>%
    mutate(
      Day = seq(0,sim$N-1,by = 1)[sample] ,
      Batch = as.factor(sim$X[2, sample]),
      Batch = ifelse(Batch == 1, "Batch1", "Batch2")
    ) %>%
    filter(!is.na(val)) %>%
    group_by(Day, coord, Batch) %>%
    summarise_posterior(val, na.rm = TRUE) %>%
    ungroup() %>%
    filter(coord == coord_keep) %>%
    dplyr::select(-coord)
}
