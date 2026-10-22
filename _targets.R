# Setup
library(targets)
library(tarchetypes)
source("R/functions_ms.R")
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
options(tidyverse.quiet = TRUE)
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# Run targets
tar_option_set(packages = c(
  "rlang", "readxl", "gtools", "rstan", "MixSIAR", "rjags",
  "R2jags", "truncnorm", "ggdist", "tools", "bayesplot", "MCMCpack",
  "R.devices", "ggsci", "ggh4x", "quarto", "knitr", "kableExtra",
  "patchwork", "cowplot", "scales", "ggridges", "viridisLite", "MixMustR",
  "tidyverse"
))
list(
  ###########################
  ###########################
  ###########################
  # MAIN TEXT
  ###########################
  ###########################
  ###########################
  tar_target(tracer_pars, MixMustR::tracer_parameters),
  tar_target(tracer_mus, tracer_pars$mus),
  tar_target(tracer_sigmas, tracer_pars$sigmas),
  tar_target(synth_df_c, MixMustR::synthetic_df_convergent),
  tar_target(synth_df_d, MixMustR::synthetic_df_divergent),
  tar_target(synth_df_c_online_si,
    write.csv(synth_df_c, "output/synthetic_dataset_1_convergent.csv")
  ),
  tar_target(synth_df_d_online_si,
    write.csv(synth_df_d, "output/synthetic_dataset_2_divergent.csv")
  ),
  tar_target(mixsiar_synth_df_c, regenarate_synth_df_iso(tracer_pars)),
  tar_target(mixsiar_objts,
    generate_mixsiar_objts(
      mixsiar_synth_df_c$pars$mus,
      mixsiar_synth_df_c$pars$sigmas,
      mixsiar_synth_df_c$df
    )
  ),
  tar_target(fig_2,
    MixMustR::compare_mixing_proportions(
      synth_df_d, synth_df_c, MixMustR::tracer_parameters$mus
    ) +
      scale_fill_manual(
        name = "Synthetic dataset:",
        labels = c("Convergent (1)", "Divergent (2)"),
        values = c("dodgerblue3", "tomato3")
      ) +
      scale_shape_manual(
        name = "Synthetic dataset:",
        labels = c("Convergent (1)", "Divergent (2)"),
        values = 21:22
      )
  ),
  tar_target(fig_3,
    plot_isospace(c(1, 2), mixsiar_objts$mix, mixsiar_objts$source,
                  mixsiar_objts$discr, mixsiar_synth_df_c$df$stream_1_props)
  ),
  tar_target(hpc_results,
    summarise_hpc_simulations(
      jobs_file = "hpc/jobs_dev.csv",
      rds_dir = "hpc/diagnostics_rds",
      out_dir = "hpc/diagnostics_rds/summaries"
    )
  ),
  tar_target(fig_4,
    make_perform_figs(
      hpc_results, truth_stream = 1, target_to_plot = "stream_1_props"
    )
  ),
  tar_target(fig_5,
    make_perform_figs(
      hpc_results, truth_stream = 2, target_to_plot = "df_stream_2"
    )
  ),
  tar_target(fig_s1,
    make_perform_figs(
      hpc_results, truth_stream = 1, target_to_plot = "df_stream_2"
    )
  ),
  tar_target(fig_s2,
    make_perform_figs(
      hpc_results, truth_stream = 2, target_to_plot = "stream_1_props"
    )
  ),
  tar_target(fig_s3,
    plot_stan_diagnostic(
      hpc_results, "n_divergent", "Divergent transitions",
      "Number of divergent transitions", trim = TRUE, limits = c(0, 2)
    )
  ),
  tar_target(fig_s4,
    plot_stan_diagnostic(
      hpc_results, "n_max_treedepth", "Max. treedepth",
      "Number of maximum-treedepth events", trim = TRUE, limits = c(0, 1)
    )
  ),
  tar_target(fig_s5,
    plot_stan_diagnostic(
      hpc_results, "min_bfmi", "Minimum BFMI", "Minimum BFMI"
    )
  ),
  tar_target(fig_s6,
    plot_stan_diagnostic(
      hpc_results, "max_rhat", "Maximum R-hat", "Maximum R-hat"
    )
  ),
  tar_target(fig_s7,
    plot_stan_diagnostic(
      hpc_results, "min_n_eff", "Min. effective sample size",
      "Minimum effective sample size"
    )
  ),
  tar_target(fig_6, prior_uncertainty_plot(x = c(0.06, 0.6, 0.04, 0.3))),
  tar_quarto(manuscript, path = "docs/MS.qmd", quiet = FALSE),
  tar_target(tbl_s1_file, "data/table_s1.csv", format = "file"),
  tar_target(tbl_s1, read.csv(tbl_s1_file, check.names = FALSE)),
  tar_quarto(supps, path = "docs/Supps.qmd", quiet = FALSE)
)
