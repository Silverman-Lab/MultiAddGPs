# Metric calculations for Simulation 1.
#
# These helpers compare model component estimates with the simulated truth.

eval_component <- function(est, truth) {
  tibble(
    RMSE = sqrt(mean((est$mean - truth$mean)^2)),
    Avg_CI_Width = mean(est$p97.5 - est$p2.5),
    Coverage = mean(truth$mean > est$p2.5 & truth$mean < est$p97.5)
  )
}

evaluate_methods <- function(estimates, truth) {
  # 1. Loop over each Model (e.g., GAM, MultiAddGPs)
  purrr::imap_dfr(estimates, function(model_est, model_name) {
    
    # 2. Get components that exist in both this model's output AND the truth
    comps <- intersect(names(model_est), names(truth))
    
    # 3. Loop over each Component (e.g., Batch, General)
    purrr::map_dfr(comps, function(comp_name) {
      
      est_df   <- model_est[[comp_name]]
      truth_df <- truth[[comp_name]]
      
      # Run the math
      eval_component(est_df, truth_df) %>%
        mutate(
          model = model_name,
          component = comp_name
        )
    })
  })
}
