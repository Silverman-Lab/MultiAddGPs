# Reproducing the MultiAddGPs analyses

R scripts for reproducing figures in the MultiAddGPs manuscript.
Use fido version 1.1.5 for all analyses and simulations.

## Dependencies

Run all commands from the repository's top-level directory. The scripts require
these R packages: `ape`, `Biostrings`, `car`, `compositions`, `dplyr`, `fido`,
`future.apply`, `ggplot2`, `ggtext`, `ggtree`, `gridExtra`, `gtable`, `lgpr`,
`MASS`, `mgcv`, `ParBayesianOptimization`, `patchwork`, `philr`, `phyloseq`,
`plotly`, `purrr`, `reshape2`, `rstan`, `scales`, `stringr`, `tibble`, and `tidyr`.
R's built-in packages, including `grid`, do not require separate installation.

## Figure filenames

Some scripts save figures under descriptive names rather than replacing the
numbered manuscript figures already supplied in `figures/`. For example:

| Generated file in `figures/` | Supplied manuscript figure |
| --- | --- |
| `final_metrics_plot.pdf` | `Main_Figure1.pdf` |
| `final_runtime_plot.pdf` | `Main_Figure2.pdf` |
| `top_balances.png` | `Main_Figure3.png` |
| `Web_Figure3_corrected.pdf` (likewise 4 and 5) | `Web_Figure3.pdf` (likewise 4 and 5) |


## Gut analyses

Human gut: Main Figure 3 and Web Figures 9–11.
Artificial gut: Main Figure 4 and Web Figures 6–8.

```sh
Rscript --vanilla Human_gut_analysis.R
Rscript --vanilla Artificial_gut_analysis.R
```
## Simulation 1

Generate Main Figures 1–2 (comparison and runtime plots) from the supplied summary tables:

```sh
Rscript --vanilla scripts/create_final_plots.R
```

To rerun the fits, prepare manifests from `config/seeds.csv`, then generate
datasets and fit models (task 1 shown):

```sh
Rscript --vanilla scripts/prepare_experiment.R
Rscript --vanilla scripts/generate_dataset_task.R manifests/datasets.csv 1
Rscript --vanilla scripts/run_model_task.R manifests/multiadd.csv 1
```
 After all fits finish:

```sh
Rscript --vanilla scripts/aggregate_experiment_results.R
Rscript --vanilla scripts/create_final_plots.R
```

## Simulation 2

Redraw Web Figures 3–5 from the saved data and fitted results:

```sh
Rscript --vanilla simulation2.R
```

To rerun hyperparameter optimization and model fitting for NAddGP and MultiAddGP:

```sh
Rscript --vanilla simulation2_fit.R
```
