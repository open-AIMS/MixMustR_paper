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
  "tidyverse", "rlang", "readxl", "gtools", "rstan", "MixSIAR", "rjags",
  "R2jags", "truncnorm", "ggdist", "tools", "bayesplot", "MCMCpack",
  "R.devices", "ggsci", "ggh4x", "quarto", "knitr", "kableExtra",
  "patchwork", "MixMustR"
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
  tar_target(mixmustr_models_c,
    MixMustR::run_mixmustr_models(
      MixMustR::mixmustr_models, synth_df_c, tracer_pars,
      sigma_ln_rho = 0.1, iter = 1e4, warmup = 5e3, chains = 4, cores = 4
    )
  ),
  tar_target(mixmustr_models_d,
    MixMustR::run_mixmustr_models(
      MixMustR::mixmustr_models, synth_df_d, tracer_pars,
      sigma_ln_rho = 1, iter = 1e4, warmup = 5e3, chains = 4, cores = 4
    )
  ),
  tar_target(mixsiar_objts,
    generate_mixsiar_objts(tracer_mus, tracer_sigmas, synth_df_c)
  ),
  tar_target(fig_2,
    plot_isospace(c(1, 2), mixsiar_objts$mix, mixsiar_objts$source,
                  mixsiar_objts$discr, synth_df_c$stream_1_props)
  ),
  tar_target(fig_3,
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
  tar_target(fig_4, {
    to_filter <- which(
      MixMustR::mixmustr_models$sample_tracer &
        MixMustR::mixmustr_models$fix_unsampled &
        MixMustR::mixmustr_models$hierarchical
    )
    a <- rbind(
      MixMustR::make_post_prop_long(
        mixmustr_models_c[[to_filter]]$model, tracer_mus, synth_df_c,
        target = "stream_1_props", n = "1 (Tracers)"
      ),
      MixMustR::make_post_prop_long(
        mixmustr_models_c[[to_filter]]$model, tracer_mus, synth_df_c,
        target = "df_stream_2", n = "2 (eDNA)"
      )
    ) |>
      MixMustR::plot_multiple_faceted_scatter_avg() +
      labs(y = "Observed (simulated)",
           title = "Model fit performance, mixing proportions",
           subtitle = "Synthetic dataset 1 (convergent)",
           fill = "Data stream:", shape = "Data stream:") +
      theme(axis.title = element_text(size = 14),
            axis.text = element_text(size = 8), 
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12), 
            legend.position = "none"
            ) +
      guides(shape = guide_legend(override.aes = list(size = 3.5)),
             fill = guide_legend(override.aes = list(alpha = 1)))
    b <- rbind(
      MixMustR::make_post_prop_long(
        mixmustr_models_d[[to_filter]]$model, tracer_mus, synth_df_d,
        target = "stream_1_props", n = "1 (Tracers)"
      ),
      MixMustR::make_post_prop_long(
        mixmustr_models_d[[to_filter]]$model, tracer_mus, synth_df_d,
        target = "df_stream_2", n = "2 (eDNA)"
      )
    ) |>
      MixMustR::plot_multiple_faceted_scatter_avg() +
      labs(y = "Observed (simulated)",
           title = "Model fit performance, mixing proportions",
           subtitle = "Synthetic dataset 2 (divergent)",
           fill = "Data stream:", shape = "Data stream:") +
      theme(axis.title = element_text(size = 14),
            axis.text = element_text(size = 8), 
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12), 
            legend.position = "inside",
            legend.position.inside = c(0.5, 0.15)
            ) +
      guides(shape = guide_legend(override.aes = list(size = 3.5)),
             fill = guide_legend(override.aes = list(alpha = 1)))
    a + b +
      patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")") +
      patchwork::plot_layout(axis_titles = "collect")
  }),
  tar_target(mixsiar_model, {
    resid_err <- TRUE
    process_err <- TRUE
    MixSIAR::write_JAGS_model(
      "MixSIAR_model_paper.txt", resid_err, process_err, mixsiar_objts$mix,
      mixsiar_objts$source
    )
    read.csv("MixSIAR_model_paper.txt")
  }),
  tar_target(jags_alpha,
    # Generate alpha hyper-parameters scaling sum(alpha)=n.sources
    synth_df_c$stream_1_props |>
      dplyr::group_by(group) |>
      dplyr::summarise(across(Macroalgae:`Terrestrial grass`, mean)) |>
      dplyr::ungroup() |>
      dplyr::select(Macroalgae:`Terrestrial grass`) |>
      colMeans() |>
      (\(x) x * length(x) / sum(x))()
  ),
  tar_target(jagsfits, {
    timing <- system.time({
      run_specs <- list(chainLength = 1e4, burn = 5e3, thin = 1, chains = 4,
                        calcDIC = TRUE)
      model <- mixsiar_run(
        run = run_specs, mixsiar_objts$mix, mixsiar_objts$source,
        mixsiar_objts$discr, "MixSIAR_model_paper.txt", alpha.prior = jags_alpha
      )
    })
    list(timing = timing, model = model)
  }),
  tar_target(tbl_1,
    make_tbl_1(mixmustr_models_c, synth_df_c, mixmustr_models_d, synth_df_d,
               jagsfits, target = "stream_1_props",
               order_ref = tracer_mus$source)
  ),
  tar_target(fig_5,
    make_jags_post_prop_long(
      jagsfits$model, dplyr::select(tracer_mus, !(`24:0`:`20:5w3`)), synth_df_c,
      target = "stream_1_props", n = 1
    ) |>
      MixMustR::plot_multiple_faceted_scatter_avg() +
      theme(legend.position = "none") +
      labs(y = "Observed (simulated) from chemical tracers (data stream 1)",
           title = "Model fit performance, mixing proportions",
           subtitle = "MixSIAR")
  ),
  tar_target(fig_6, prior_uncertainty_plot(x = c(0.06, 0.6, 0.04, 0.3))),
  tar_quarto(manuscript, path = "docs/MS.qmd", quiet = FALSE),
  ###########################
  ###########################
  ###########################
  # SUPPLEMENTARY INFORMATION
  # AND DIAGNOSTICS
  ###########################
  ###########################
  ###########################
  tar_target(fig_s1,
    make_fig_s1_8(
      to_filter = 2:4, model_list = mixmustr_models_c,
      ylab = "Observed (simulated) from chemical tracers (data stream 1)",
      subt = "Synthetic dataset 1 (convergent)",
      mu_tab = tracer_mus, data_streams_list = synth_df_c,
      target = "stream_1_props"
    )
  ),
  tar_target(fig_s2,
    make_fig_s1_8(
      to_filter = 5:8, model_list = mixmustr_models_c,
      ylab = "Observed (simulated) from chemical tracers (data stream 1)",
      subt = "Synthetic dataset 1 (convergent)",
      mu_tab = tracer_mus, data_streams_list = synth_df_c,
      target = "stream_1_props")
  ),
  tar_target(fig_s3,
    make_fig_s1_8(
      to_filter = 2:4, model_list = mixmustr_models_c,
      ylab = "Observed (simulated) from eDNA (data stream 2)",
      subt = "Synthetic dataset 1 (convergent)",
      mu_tab = tracer_mus, data_streams_list = synth_df_c,
      target = "df_stream_2"
    )
  ),
  tar_target(fig_s4,
    make_fig_s1_8(
      to_filter = 5:8, model_list = mixmustr_models_c,
      ylab = "Observed (simulated) from eDNA (data stream 2)",
      subt = "Synthetic dataset 1 (convergent)",
      mu_tab = tracer_mus, data_streams_list = synth_df_c,
      target = "df_stream_2"
    )
  ),
  tar_target(fig_s5,
    make_fig_s1_8(
      to_filter = 2:4, model_list = mixmustr_models_d,
      ylab = "Observed (simulated) from chemical tracers (data stream 1)",
      subt = "Synthetic dataset 2 (divergent)",
      mu_tab = tracer_mus, data_streams_list = synth_df_d,
      target = "stream_1_props"
    )
  ),
  tar_target(fig_s6,
    make_fig_s1_8(
      to_filter = 5:8, model_list = mixmustr_models_d,
      ylab = "Observed (simulated) from chemical tracers (data stream 1)",
      subt = "Synthetic dataset 2 (divergent)",
      mu_tab = tracer_mus, data_streams_list = synth_df_d,
      target = "stream_1_props")
  ),
  tar_target(fig_s7,
    make_fig_s1_8(
      to_filter = 2:4, model_list = mixmustr_models_d,
      ylab = "Observed (simulated) from eDNA (data stream 2)",
      subt = "Synthetic dataset 2 (divergent)",
      mu_tab = tracer_mus, data_streams_list = synth_df_d,
      target = "df_stream_2"
    )
  ),
  tar_target(fig_s8,
    make_fig_s1_8(
      to_filter = 5:8, model_list = mixmustr_models_d,
      ylab = "Observed (simulated) from eDNA (data stream 2)",
      subt = "Synthetic dataset 2 (divergent)",
      mu_tab = tracer_mus, data_streams_list = synth_df_d,
      target = "df_stream_2"
    )
  ),
  tar_target(fig_s9, plot_supp_checks(mixmustr_models_c[[1]]$model, 1)),
  tar_target(fig_s10, plot_supp_checks(mixmustr_models_c[[2]]$model, 2)),
  tar_target(fig_s11, plot_supp_checks(mixmustr_models_c[[3]]$model, 3)),
  tar_target(fig_s12, plot_supp_checks(mixmustr_models_c[[4]]$model, 4)),
  tar_target(fig_s13, plot_supp_checks(mixmustr_models_c[[5]]$model, 5)),
  tar_target(fig_s14, plot_supp_checks(mixmustr_models_c[[6]]$model, 6)),
  tar_target(fig_s15, plot_supp_checks(mixmustr_models_c[[7]]$model, 7)),
  tar_target(fig_s16, plot_supp_checks(mixmustr_models_c[[8]]$model, 8)),
  tar_target(fig_s17, plot_supp_checks(mixmustr_models_d[[1]]$model, 1)),
  tar_target(fig_s18, plot_supp_checks(mixmustr_models_d[[2]]$model, 2)),
  tar_target(fig_s19, plot_supp_checks(mixmustr_models_d[[3]]$model, 3)),
  tar_target(fig_s20, plot_supp_checks(mixmustr_models_d[[4]]$model, 4)),
  tar_target(fig_s21, plot_supp_checks(mixmustr_models_d[[5]]$model, 5)),
  tar_target(fig_s22, plot_supp_checks(mixmustr_models_d[[6]]$model, 6)),
  tar_target(fig_s23, plot_supp_checks(mixmustr_models_d[[7]]$model, 7)),
  tar_target(fig_s24, plot_supp_checks(mixmustr_models_d[[8]]$model, 8)),
  tar_target(tbl_s1_file, "data/table_s1.csv", format = "file"),
  tar_target(tbl_s1, read.csv(tbl_s1_file, check.names = FALSE)),
  tar_target(tbl_s2, make_table_s2(mixmustr_models_c, mixmustr_models_d)),

  tar_target(prior_model_path,
    "stan/sample_tracer_FALSE_fix_unsampled_TRUE_hierarchical_FALSE_ln_rho_as_prior.stan", format = "file"
  ),
  tar_target(prior_model,
    run_prior_model_si(synth_df_d, tracer_pars, prior_model_path)
  ),
  tar_target(fig_s26,
    make_fig_s1_8(
      to_filter = 6, model_list = prior_model,
      ylab = "Observed (simulated) from chemical tracers (data stream 1)",
      subt = "Synthetic dataset 2 (divergent)",
      mu_tab = tracer_mus, data_streams_list = synth_df_d,
      target = "stream_1_props"
    )
  ),
  tar_target(fig_s27,
    make_fig_s1_8(
      to_filter = 6, model_list = prior_model,
      ylab = "Observed (simulated) from eDNA (data stream 2)",
      subt = "Synthetic dataset 2 (divergent)",
      mu_tab = tracer_mus, data_streams_list = synth_df_d,
      target = "df_stream_2"
    )
  ),
  tar_target(fig_s28, plot_supp_checks(prior_model[[6]]$model, 6)),
  tar_quarto(supps, path = "docs/Supps.qmd", quiet = FALSE)
)
