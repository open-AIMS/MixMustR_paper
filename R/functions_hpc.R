# MixMustR paper per-simulation diagnostics pipeline
# -----------------------------------------------------------------------------
# This pipeline is designed for the MixMustR workflow where each
# simulation fits one row of `MixMustR::mixmustr_models` via:
#
#   model_run <- MixMustR::run_mixmustr_models(
#     MixMustR::mixmustr_models[i, ], synth_df, tracer_parameters, bcs,
#     sigma_ln_rho = 1, iter = 2e3, warmup = 1e3,
#     chains = 4, cores = 4, seed = stan_seed,
#     control = list(adapt_delta = 0.99, max_treedepth = 15)
#   )
#
# The key design goal is to compute lightweight diagnostics immediately after
# each model fit, save those diagnostics, and save the heavy Stan model
# elsewhere. At the end of all simulations, the saved lightweight diagnostics
# can be aggregated across simulations and model variants.
#
# Suggested simulation-loop usage:
#
#   model_choice <- MixMustR::mixmustr_models[model_id, ]
#
#   model_run <- MixMustR::run_mixmustr_models(
#     model_choice, synth_df, tracer_parameters, bcs,
#     sigma_ln_rho = 1, iter = 2e3, warmup = 1e3,
#     chains = 4, cores = 4, seed = stan_seed
#   )
#
#   diag <- diagnose_mixmustr_model_run(
#     model_run = model_run,
#     model_choice = model_choice,
#     synth_df = synth_df,
#     tracer_list = tracer_parameters,
#     sign_df = bcs,
#     sigma_ln_rho = 1,
#     sim_id = sim_id,
#     prop_seed = prop_seed,
#     tracer_seed = tracer_seed,
#     stan_seed = stan_seed,
#     output_dir = "diagnostics_rds",
#     max_treedepth = 15L
#   )
#
#   rm(model_run); gc()
#
# End-of-simulation aggregation:
#
#   diagnostics <- aggregate_mixmustr_diagnostics(
#     diag_files = list.files(
#       "diagnostics_rds", pattern = "\\.rds$", full.names = TRUE
#     ), output_dir = "diagnostics_csv",
#     write_csv = TRUE
#   )
#   variant_summary <- summarise_mixmustr_diagnostics_by_variant(diagnostics)

`%||%` <- function(x, y) if (is.null(x)) y else x

is_stanfit <- function(x) inherits(x, "stanfit")

build_variant_label <- function(model_choice, model_id = NA_integer_) {
  st <- as.logical(model_choice$sample_tracer[[1]])
  fu <- as.logical(model_choice$fix_unsampled[[1]])
  hi <- as.logical(model_choice$hierarchical[[1]])

  paste0(
    "model_", ifelse(is.na(model_id), NA_character_, model_id),
    "__sample_tracer_", toupper(as.character(st)),
    "__fix_unsampled_", toupper(as.character(fu)),
    "__hierarchical_", toupper(as.character(hi))
  )
}

extract_fit_from_model_run <- function(model_run) {
  if (is_stanfit(model_run)) return(model_run)

  if (is.list(model_run) && length(model_run) >= 1) {
    if ("model" %in% names(model_run) && is_stanfit(model_run$model)) return(model_run$model)
    if (is.list(model_run[[1]]) && "model" %in% names(model_run[[1]]) && is_stanfit(model_run[[1]]$model)) {
      return(model_run[[1]]$model)
    }
  }

  stop("Could not locate a stanfit object in `model_run`.", call. = FALSE)
}

stan_param_exists <- function(fit, par) {
  par %in% fit@sim$pars_oi
}

bfmi_one_chain <- function(energy) {
  energy <- as.numeric(energy)
  if (length(energy) < 3 || stats::var(energy) == 0) return(NA_real_)
  mean(diff(energy)^2) / stats::var(energy)
}

# Exclude deterministic sampled-source rows of all_mus from scalar Stan
# diagnostics when sample_tracer = FALSE.
#
# Why:
#   In variants with sample_tracer = FALSE, all_mus[1:J, k] are deterministic
#   transformed parameters equal to x[1:J, k]. rstan::summary() still reports
#   n_eff/Rhat for them, and those diagnostics can dominate min_n_eff even
#   though no sampling is occurring for those quantities.
#
# Policy:
#   * If sample_tracer = FALSE: drop all_mus[j, k] for j = 1, ..., J from
#     max_rhat/min_n_eff summaries.
#   * If fix_unsampled = FALSE: keep all_mus[J + 1, k], because it corresponds
#     to the estimated unsampled signature.
#   * Optional: if exclude_fixed_unsampled = TRUE and fix_unsampled = TRUE, also
#     drop all_mus[J + 1, k], because that row is deterministic too.

params_for_diagnostics <- function(param_names,
                                   sdata = NULL,
                                   sample_tracer = TRUE,
                                   fix_unsampled = FALSE,
                                   exclude_fixed_unsampled = TRUE) {
  keep <- rep(TRUE, length(param_names))

  if (!isTRUE(sample_tracer)) {
    J <- NULL

    if (!is.null(sdata) && "J" %in% names(sdata)) {
      J <- as.integer(sdata$J)
    }

    # Fallback: infer J as max first index in all_mus minus one, assuming
    # all_mus has dimensions (J + 1) x K.
    if (is.null(J) || is.na(J)) {
      idx <- regmatches(
        param_names,
        regexec("^all_mus\\[([0-9]+),([0-9]+)\\]$", param_names)
      )
      first_idx <- suppressWarnings(as.integer(vapply(
        idx,
        function(z) if (length(z) == 3) z[[2]] else NA_character_,
        character(1)
      )))

      if (any(is.finite(first_idx))) {
        J <- max(first_idx, na.rm = TRUE) - 1L
      }
    }

    if (!is.null(J) && is.finite(J) && J >= 1L) {
      idx <- regmatches(
        param_names,
        regexec("^all_mus\\[([0-9]+),([0-9]+)\\]$", param_names)
      )
      first_idx <- suppressWarnings(as.integer(vapply(
        idx,
        function(z) if (length(z) == 3) z[[2]] else NA_character_,
        character(1)
      )))

      sampled_all_mus <- !is.na(first_idx) & first_idx <= J
      keep[sampled_all_mus] <- FALSE
    }
  }

  if (isTRUE(exclude_fixed_unsampled) && isTRUE(fix_unsampled)) {
    J <- if (!is.null(sdata) && "J" %in% names(sdata)) as.integer(sdata$J) else NA_integer_
    if (is.finite(J)) {
      unsampled_pattern <- paste0("^all_mus\\[", J + 1L, ",")
      keep[grepl(unsampled_pattern, param_names)] <- FALSE
    }
  }

  keep
}

summarise_stan_diagnostics <- function(fit,
                                       metadata,
                                       sdata = NULL,
                                       max_treedepth = 15L,
                                       exclude_fixed_unsampled = TRUE) {
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)

  n_divergent <- sum(vapply(
    sp,
    function(x) sum(x[, "divergent__"]),
    numeric(1)
  ))

  n_max_treedepth <- sum(vapply(
    sp,
    function(x) sum(x[, "treedepth__"] >= max_treedepth),
    numeric(1)
  ))

  bfmi_one_chain <- function(energy) {
    energy <- as.numeric(energy)
    if (length(energy) < 3 || stats::var(energy) == 0) return(NA_real_)
    mean(diff(energy)^2) / stats::var(energy)
  }

  bfmi <- vapply(sp, function(x) bfmi_one_chain(x[, "energy__"]), numeric(1))
  min_bfmi <- suppressWarnings(min(bfmi, na.rm = TRUE))
  if (!is.finite(min_bfmi)) min_bfmi <- NA_real_

  summ <- rstan::summary(fit)$summary

  keep_param <- params_for_diagnostics(
    param_names = rownames(summ),
    sdata = sdata,
    sample_tracer = metadata$sample_tracer[[1]],
    fix_unsampled = metadata$fix_unsampled[[1]],
    exclude_fixed_unsampled = exclude_fixed_unsampled
  )

  summ_diag <- summ[keep_param, , drop = FALSE]
  rhat <- summ_diag[, "Rhat"]
  n_eff <- summ_diag[, "n_eff"]

  dplyr::bind_cols(
    metadata,
    tibble::tibble(
      n_divergent = n_divergent,
      n_max_treedepth = n_max_treedepth,
      min_bfmi = min_bfmi,
      max_rhat = suppressWarnings(max(rhat[is.finite(rhat)], na.rm = TRUE)),
      min_n_eff = suppressWarnings(min(n_eff[is.finite(n_eff)], na.rm = TRUE)),
      n_params_reported = nrow(summ),
      n_params_in_diagnostics = nrow(summ_diag),
      n_params_excluded_from_diagnostics = sum(!keep_param)
    )
  )
}

reconstruct_sdata_from_dev_wrangle <- function(synth_df,
                                               tracer_list,
                                               model_choice,
                                               sigma_ln_rho,
                                               sign_df = NULL) {
  wrangle_fun <- get("mixmustr_wrangle_input", envir = asNamespace("MixMustR"))

  args <- list(
    data_streams_list = synth_df,
    tracer_list = tracer_list,
    model_path = model_choice$code_path[[1]] %||% NULL,
    sigma_ln_rho = sigma_ln_rho,
    sample_tracer = model_choice$sample_tracer[[1]],
    fix_unsampled = model_choice$fix_unsampled[[1]],
    hierarchical = model_choice$hierarchical[[1]]
  )

  # Keep this flexible in case the dev branch adds/removes optional arguments.
  fmls <- names(formals(wrangle_fun))
  args <- args[names(args) %in% fmls]

  do.call(wrangle_fun, args)
}

summarise_prop_recovery <- function(post_long,
                                    source_col = "source",
                                    expected_col = "Predicted",
                                    observed_col = "Observed") {
  if (is.null(post_long) || !all(c(source_col, expected_col, observed_col) %in% names(post_long))) {
    return(tibble::tibble())
  }

  post_long |>
    dplyr::rename(
      source = dplyr::all_of(source_col),
      predicted = dplyr::all_of(expected_col),
      observed = dplyr::all_of(observed_col)
    ) |>
    dplyr::mutate(
      error = .data$predicted - .data$observed,
      abs_error = abs(.data$error),
      sq_error = .data$error^2
    ) |>
    dplyr::group_by(.data$source) |>
    dplyr::summarise(
      n = dplyr::n(),
      rmse = sqrt(mean(.data$sq_error, na.rm = TRUE)),
      mae = mean(.data$abs_error, na.rm = TRUE),
      bias = mean(.data$error, na.rm = TRUE),
      cor = suppressWarnings(stats::cor(.data$predicted, .data$observed, use = "complete.obs")),
      .groups = "drop"
    )
}

summarise_prop_recovery_target <- function(fit, tracer_list, synth_df, target, metadata) {
  post_long <- tryCatch(
    MixMustR::make_post_prop_long(
      fit,
      tracer_list$mus,
      synth_df,
      target = target,
      n = 1
    ),
    error = function(e) NULL
  )

  out <- summarise_prop_recovery(post_long)
  if (!nrow(out)) return(out)

  dplyr::bind_cols(metadata, tibble::tibble(target = target), out)
}

get_tracer_center_scale <- function(sdata) {
  Y <- as.matrix(sdata$Y)
  list(
    center = attr(sdata, "tracer_center") %||% rep(0, ncol(Y)),
    scale = attr(sdata, "tracer_scale") %||% rep(1, ncol(Y))
  )
}

backtransform_draw_matrix <- function(draw_mat, sdata) {
  cs <- get_tracer_center_scale(sdata)
  out <- sweep(draw_mat, 2, cs$scale, "*")
  sweep(out, 2, cs$center, "+")
}

summarise_true_residuals <- function(synth_df, tracer_list, sign_df, metadata) {
  if (is.null(sign_df)) return(tibble::tibble())
  if (!all(c("df_stream_1", "stream_1_props") %in% names(synth_df))) return(tibble::tibble())
  if (!"Unsampled" %in% sign_df$source) return(tibble::tibble())

  tracer_names <- setdiff(names(tracer_list$mus), "source")

  S_sampled <- tracer_list$mus |>
    tibble::column_to_rownames("source") |>
    as.matrix()

  S_unsampled <- sign_df |>
    dplyr::filter(.data$source == "Unsampled") |>
    dplyr::select(dplyr::all_of(colnames(S_sampled))) |>
    as.matrix()

  S_all <- rbind(S_sampled, Unsampled = S_unsampled[1, ])

  P_sampled <- synth_df$stream_1_props |>
    dplyr::select(dplyr::all_of(tracer_list$mus$source)) |>
    as.matrix()

  P_all <- cbind(P_sampled, Unsampled = 1 - rowSums(P_sampled))
  Y_hat <- P_all %*% S_all

  Y_obs <- synth_df$df_stream_1 |>
    dplyr::select(dplyr::all_of(tracer_names)) |>
    as.matrix()

  resid_true <- Y_obs - Y_hat

  dplyr::bind_cols(
    metadata,
    tibble::tibble(
      tracer = colnames(resid_true),
      true_resid_mean = as.numeric(colMeans(resid_true, na.rm = TRUE)),
      true_resid_sd = as.numeric(apply(resid_true, 2, stats::sd, na.rm = TRUE)),
      true_resid_min = as.numeric(apply(resid_true, 2, min, na.rm = TRUE)),
      true_resid_max = as.numeric(apply(resid_true, 2, max, na.rm = TRUE))
    )
  )
}

summarise_total_sigma <- function(fit, sdata, metadata) {
  if (is.null(sdata) || !all(c("s", "Y") %in% names(sdata))) return(tibble::tibble())
  if (!stan_param_exists(fit, "p") || !stan_param_exists(fit, "sigma_y")) return(tibble::tibble())

  draws <- rstan::extract(fit, pars = c("p", "sigma_y"), permuted = TRUE)
  p_draws <- draws$p
  sigma_y <- draws$sigma_y
  s <- as.matrix(sdata$s)

  n_iter <- dim(p_draws)[1]
  N <- dim(p_draws)[2]
  J <- nrow(s)
  K <- ncol(s)

  total_sigma_std <- matrix(NA_real_, nrow = n_iter, ncol = K)

  for (d in seq_len(n_iter)) {
    for (k in seq_len(K)) {
      process_var_n <- rep(sigma_y[d, k]^2, N)
      for (j in seq_len(J)) {
        process_var_n <- process_var_n + p_draws[d, , j]^2 * s[j, k]^2
      }
      total_sigma_std[d, k] <- sqrt(mean(process_var_n))
    }
  }

  tracer_scale <- get_tracer_center_scale(sdata)$scale
  total_sigma_orig <- sweep(total_sigma_std, 2, tracer_scale, "*")

  tracer_names <- colnames(as.matrix(sdata$Y)) %||% paste0("tracer_", seq_len(ncol(total_sigma_orig)))
  colnames(total_sigma_orig) <- tracer_names

  out <- as.data.frame(t(apply(
    total_sigma_orig,
    2,
    stats::quantile,
    probs = c(0.025, 0.5, 0.975),
    na.rm = TRUE
  ))) |>
    tibble::rownames_to_column("tracer") |>
    dplyr::rename(q2.5 = `2.5%`, median = `50%`, q97.5 = `97.5%`)

  dplyr::bind_cols(metadata, out)
}

summarise_mu_unsampled <- function(fit, sdata, sign_df, metadata) {
  if (!stan_param_exists(fit, "mu_unsampled")) return(tibble::tibble())
  if (is.null(sdata)) return(tibble::tibble())

  draws <- rstan::extract(fit, pars = "mu_unsampled", permuted = TRUE)$mu_unsampled
  draws_orig <- backtransform_draw_matrix(draws, sdata)

  tracer_names <- colnames(as.matrix(sdata$Y)) %||% paste0("tracer_", seq_len(ncol(draws_orig)))
  colnames(draws_orig) <- tracer_names

  out <- as.data.frame(t(apply(
    draws_orig,
    2,
    stats::quantile,
    probs = c(0.025, 0.5, 0.975),
    na.rm = TRUE
  ))) |>
    tibble::rownames_to_column("tracer") |>
    dplyr::rename(q2.5 = `2.5%`, median = `50%`, q97.5 = `97.5%`)

  if (!is.null(sign_df) && "source" %in% names(sign_df) && "Unsampled" %in% sign_df$source) {
    true_u <- sign_df |>
      dplyr::filter(.data$source == "Unsampled") |>
      dplyr::select(dplyr::all_of(tracer_names)) |>
      tidyr::pivot_longer(dplyr::everything(), names_to = "tracer", values_to = "truth")

    out <- out |>
      dplyr::left_join(true_u, by = "tracer") |>
      dplyr::mutate(
        covered = .data$truth >= .data$q2.5 & .data$truth <= .data$q97.5,
        abs_error_median = abs(.data$median - .data$truth)
      )
  }

  dplyr::bind_cols(metadata, out)
}

summarise_mu_source <- function(fit, sdata, tracer_list, metadata) {
  if (!stan_param_exists(fit, "mu_source")) return(tibble::tibble())
  if (is.null(sdata)) return(tibble::tibble())

  draws <- rstan::extract(fit, pars = "mu_source", permuted = TRUE)$mu_source
  cs <- get_tracer_center_scale(sdata)

  tracer_names <- colnames(as.matrix(sdata$Y)) %||% paste0("tracer_", seq_len(dim(draws)[3]))
  source_names <- if (!is.null(tracer_list)) {
    as.character(tracer_list$mus$source)
  } else {
    paste0("source_", seq_len(dim(draws)[2]))
  }

  out <- purrr::map_dfr(seq_along(source_names), function(j) {
    mat_std <- draws[, j, , drop = FALSE][, 1, ]
    mat_orig <- sweep(mat_std, 2, cs$scale, "*") |>
      sweep(2, cs$center, "+")
    colnames(mat_orig) <- tracer_names

    as.data.frame(t(apply(
      mat_orig,
      2,
      stats::quantile,
      probs = c(0.025, 0.5, 0.975),
      na.rm = TRUE
    ))) |>
      tibble::rownames_to_column("tracer") |>
      dplyr::rename(q2.5 = `2.5%`, median = `50%`, q97.5 = `97.5%`) |>
      dplyr::mutate(source = source_names[j], .before = 1)
  })

  if (!is.null(tracer_list)) {
    truth <- tracer_list$mus |>
      tidyr::pivot_longer(-dplyr::all_of("source"), names_to = "tracer", values_to = "truth")

    out <- out |>
      dplyr::left_join(truth, by = c("source", "tracer")) |>
      dplyr::mutate(
        covered = .data$truth >= .data$q2.5 & .data$truth <= .data$q97.5,
        abs_error_median = abs(.data$median - .data$truth)
      )
  }

  dplyr::bind_cols(metadata, out)
}

model_metadata_tbl <- function(model_choice,
                               sim_id = NA_integer_,
                               model_id = NA_integer_,
                               prop_seed = NA_integer_,
                               tracer_seed = NA_integer_,
                               stan_seed = NA_integer_) {
  variant <- build_variant_label(model_choice, model_id = model_id)

  tibble::tibble(
    sim_id = sim_id,
    model_id = model_id,
    variant = variant,
    sample_tracer = as.logical(model_choice$sample_tracer[[1]]),
    fix_unsampled = as.logical(model_choice$fix_unsampled[[1]]),
    hierarchical = as.logical(model_choice$hierarchical[[1]]),
    prop_seed = prop_seed,
    tracer_seed = tracer_seed,
    stan_seed = stan_seed
  )
}

diagnose_mixmustr_model_run <- function(model_run,
                                        model_choice,
                                        synth_df,
                                        tracer_list,
                                        sign_df = NULL,
                                        sigma_ln_rho = 1,
                                        sim_id = NA_integer_,
                                        model_id = NA_integer_,
                                        prop_seed = NA_integer_,
                                        tracer_seed = NA_integer_,
                                        stan_seed = NA_integer_,
                                        max_treedepth = 15L,
                                        output_dir = NULL,
                                        save_rds = !is.null(output_dir)) {
  fit <- extract_fit_from_model_run(model_run)

  metadata <- model_metadata_tbl(
    model_choice = model_choice,
    sim_id = sim_id,
    model_id = model_id,
    prop_seed = prop_seed,
    tracer_seed = tracer_seed,
    stan_seed = stan_seed
  )

  sdata <- if (is.list(model_run) && "data" %in% names(model_run) && !is.null(model_run$data)) {
    model_run$data
  } else {
    tryCatch(
      reconstruct_sdata_from_dev_wrangle(
        synth_df = synth_df, tracer_list = tracer_list,
        model_choice = model_choice, sigma_ln_rho = sigma_ln_rho,
        sign_df = sign_df
      ), error = function(e) NULL
    )
  }
  stan_diag <- summarise_stan_diagnostics(
    fit = fit,
    metadata = metadata,
    sdata = sdata,
    max_treedepth = max_treedepth,
    exclude_fixed_unsampled = TRUE
  )

  prop_stream1 <- summarise_prop_recovery_target(
    fit = fit,
    tracer_list = tracer_list,
    synth_df = synth_df,
    target = "stream_1_props",
    metadata = metadata
  )

  prop_stream2 <- summarise_prop_recovery_target(
    fit = fit,
    tracer_list = tracer_list,
    synth_df = synth_df,
    target = "df_stream_2",
    metadata = metadata
  )

  true_resid <- summarise_true_residuals(
    synth_df = synth_df,
    tracer_list = tracer_list,
    sign_df = sign_df,
    metadata = metadata
  )

  sigma_total <- summarise_total_sigma(fit, sdata, metadata)

  if (nrow(sigma_total) && nrow(true_resid)) {
    sigma_total <- sigma_total |>
      dplyr::left_join(
        true_resid |>
          dplyr::select("sim_id", "model_id", "variant", "tracer", "true_resid_sd"),
        by = c("sim_id", "model_id", "variant", "tracer")
      ) |>
      dplyr::mutate(
        ratio = .data$median / .data$true_resid_sd,
        abs_log_ratio = abs(log(.data$ratio))
      )
  }

  out <- list(
    stan_diagnostics = stan_diag,
    prop_recovery = dplyr::bind_rows(prop_stream1, prop_stream2),
    true_resid_summary = true_resid,
    total_sigma = sigma_total,
    mu_unsampled = summarise_mu_unsampled(fit, sdata, sign_df, metadata),
    mu_source = summarise_mu_source(fit, sdata, tracer_list, metadata)
  )

  if (isTRUE(save_rds)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    file_stub <- paste0(
      "diag_model_", metadata$model_id[[1]],
      "_sim_", metadata$sim_id[[1]],
      ".rds"
    )
    saveRDS(out, file.path(output_dir, file_stub))
  }

  out
}

safe_diagnose_mixmustr_model_run <- function(...) {
  tryCatch(
    diagnose_mixmustr_model_run(...),
    error = function(e) {
      dots <- list(...)
      model_choice <- dots$model_choice
      metadata <- model_metadata_tbl(
        model_choice = model_choice,
        sim_id = dots$sim_id %||% NA_integer_,
        model_id = dots$model_id %||% NA_integer_,
        prop_seed = dots$prop_seed %||% NA_integer_,
        tracer_seed = dots$tracer_seed %||% NA_integer_,
        stan_seed = dots$stan_seed %||% NA_integer_
      )

      list(
        stan_diagnostics = dplyr::bind_cols(
          metadata,
          tibble::tibble(
            n_divergent = NA_integer_,
            n_max_treedepth = NA_integer_,
            min_bfmi = NA_real_,
            max_rhat = NA_real_,
            min_n_eff = NA_real_,
            error_message = conditionMessage(e)
          )
        ),
        prop_recovery = tibble::tibble(),
        true_resid_summary = tibble::tibble(),
        total_sigma = tibble::tibble(),
        mu_unsampled = tibble::tibble(),
        mu_source = tibble::tibble()
      )
    }
  )
}

aggregate_mixmustr_diagnostics <- function(diag_files = NULL,
                                           diag_list = NULL,
                                           output_dir = NULL,
                                           write_csv = FALSE) {
  if (!is.null(diag_files)) {
    diag_list <- lapply(diag_files, readRDS)
  }

  if (is.null(diag_list)) {
    stop("Provide either `diag_files` or `diag_list`.", call. = FALSE)
  }

  out <- list(
    stan_diagnostics = purrr::map_dfr(diag_list, "stan_diagnostics"),
    prop_recovery = purrr::map_dfr(diag_list, "prop_recovery"),
    true_resid_summary = purrr::map_dfr(diag_list, "true_resid_summary"),
    total_sigma = purrr::map_dfr(diag_list, "total_sigma"),
    mu_unsampled = purrr::map_dfr(diag_list, "mu_unsampled"),
    mu_source = purrr::map_dfr(diag_list, "mu_source")
  )

  if (isTRUE(write_csv)) {
    if (is.null(output_dir)) stop("`output_dir` is required when `write_csv = TRUE`.", call. = FALSE)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    purrr::iwalk(out, function(x, nm) {
      readr::write_csv(x, file.path(output_dir, paste0("mixmustr_", nm, ".csv")))
    })
  }

  out
}

summarise_mixmustr_diagnostics_by_variant <- function(diagnostics) {
  list(
    stan = diagnostics$stan_diagnostics |>
      dplyr::group_by(.data$model_id, .data$variant, .data$sample_tracer, .data$fix_unsampled, .data$hierarchical) |>
      dplyr::summarise(
        n_sims = dplyr::n(),
        total_divergences = sum(.data$n_divergent, na.rm = TRUE),
        sims_with_divergences = sum(.data$n_divergent > 0, na.rm = TRUE),
        total_max_treedepth = sum(.data$n_max_treedepth, na.rm = TRUE),
        sims_with_max_treedepth = sum(.data$n_max_treedepth > 0, na.rm = TRUE),
        min_bfmi = min(.data$min_bfmi, na.rm = TRUE),
        max_rhat = max(.data$max_rhat, na.rm = TRUE),
        sims_rhat_gt_1.01 = sum(.data$max_rhat > 1.01, na.rm = TRUE),
        .groups = "drop"
      ),

    prop_recovery = diagnostics$prop_recovery |>
      dplyr::group_by(.data$model_id, .data$variant, .data$target, .data$source) |>
      dplyr::summarise(
        median_rmse = stats::median(.data$rmse, na.rm = TRUE),
        mean_rmse = mean(.data$rmse, na.rm = TRUE),
        q95_rmse = stats::quantile(.data$rmse, 0.95, na.rm = TRUE),
        median_mae = stats::median(.data$mae, na.rm = TRUE),
        mean_abs_bias = mean(abs(.data$bias), na.rm = TRUE),
        median_cor = stats::median(.data$cor, na.rm = TRUE),
        min_cor = min(.data$cor, na.rm = TRUE),
        .groups = "drop"
      ),

    total_sigma = diagnostics$total_sigma |>
      dplyr::group_by(.data$model_id, .data$variant, .data$tracer) |>
      dplyr::summarise(
        median_ratio = stats::median(.data$ratio, na.rm = TRUE),
        mean_ratio = mean(.data$ratio, na.rm = TRUE),
        q05_ratio = stats::quantile(.data$ratio, 0.05, na.rm = TRUE),
        q95_ratio = stats::quantile(.data$ratio, 0.95, na.rm = TRUE),
        median_abs_log_ratio = stats::median(.data$abs_log_ratio, na.rm = TRUE),
        .groups = "drop"
      ),

    mu_unsampled = diagnostics$mu_unsampled |>
      dplyr::group_by(.data$model_id, .data$variant, .data$tracer) |>
      dplyr::summarise(
        coverage = mean(.data$covered, na.rm = TRUE),
        median_abs_error = stats::median(.data$abs_error_median, na.rm = TRUE),
        mean_abs_error = mean(.data$abs_error_median, na.rm = TRUE),
        q95_abs_error = stats::quantile(.data$abs_error_median, 0.95, na.rm = TRUE),
        .groups = "drop"
      ),

    mu_source = diagnostics$mu_source |>
      dplyr::group_by(.data$model_id, .data$variant, .data$source, .data$tracer) |>
      dplyr::summarise(
        coverage = mean(.data$covered, na.rm = TRUE),
        median_abs_error = stats::median(.data$abs_error_median, na.rm = TRUE),
        mean_abs_error = mean(.data$abs_error_median, na.rm = TRUE),
        q95_abs_error = stats::quantile(.data$abs_error_median, 0.95, na.rm = TRUE),
        .groups = "drop"
      )
  )
}

simulate_mvn_mixture <- function(G, M, J, noise = 0.5, seed = 1) {
  N <- G * M
  # 1. Softmax function to convert real numbers to proportions
  softmax <- function(x) {
    exp_x <- exp(x)
    return(exp_x / sum(exp_x))
  }
  # 2. Setup Group Centers and Covariance
  # Each group has a unique mean vector in J-dimensional space
  set.seed(seed)
  group_means_unconstrained <- matrix(
    rnorm(G * J, mean = 0, sd = 2), nrow = G, ncol = J
  )
  # Shared covariance matrix (Identity * noise)
  # You could customize this to add correlations between specific sources
  sigma <- diag(J) * noise
  # 3. Generate Data
  results_list <- list()
  for(g in 1:G) {
    # Generate M observations for group g in unconstrained space
    unconstrained_data <- MASS::mvrnorm(
      n = M, mu = group_means_unconstrained[g, ], Sigma = sigma
    )
    # Map to simplex (proportions) row by row
    proportions_data <- t(apply(unconstrained_data, 1, softmax))
    # Store with group ID
    group_df <- data.frame(group = rep(g, M))
    group_df <- cbind(group_df, proportions_data)
    results_list[[g]] <- group_df
  }
  # 4. Combine and Label
  full_dataset <- do.call(rbind, results_list)
  colnames(full_dataset)[2:(J+1)] <- paste0("Source_", 1:J)
  full_dataset$group <- factor(full_dataset$group)
  # 5. Calculate Group-Source Means
  group_means_table <- aggregate(. ~ group, data = full_dataset, FUN = mean)
  list(dataset = full_dataset, group_means = group_means_table)
}
make_mixture_data_si <- function(si_df, stream_1_props, stream_2_props,
                                 truth_stream = 1, template_df, ...) {
  if (truth_stream == 1) {
    mixing_props_df <- as.matrix(stream_1_props)
  } else if (truth_stream == 2) {
    mixing_props_df <- as.matrix(stream_2_props)
  } else {
    stop("truth_stream can only have value 1 or 2.")
  }
  out <- MixMustR:::reshape_isotope_df(si_df) |>
    dplyr::mutate(
      a = dplyr::recode_values(tracer, "d(13C/12C)" ~ -Inf, default = 0),
      b = Inf
    ) |>
    MixMustR:::calc_tracer_estimate(...) |>
    dplyr::arrange(source, tracer_family, tracer) |> # ensure alphabetical order
    split(f = ~ tracer_family + tracer, drop = TRUE) |>
    purrr::map(function(tracer_df, props_df) {
      (t(props_df[, tracer_df$source]) * tracer_df$estimate) |>
        colSums() |>
        data.frame(check.names = FALSE) |>
        MixMustR:::assign_new_names(unique(tracer_df$tracer))
    }, props_df = mixing_props_df)
  # check if bind_cols is working as intended
  df_stream_1 <- cbind(template_df, dplyr::bind_cols(out))
  df_stream_2 <- cbind(template_df, stream_2_props)
  stream_1_props <- cbind(template_df, stream_1_props)
  list(df_stream_1 = df_stream_1, df_stream_2 = df_stream_2,
       stream_1_props = stream_1_props)
}
rm_unsampled <- function(x) {
  x$df_stream_2 <- dplyr::select(x$df_stream_2, !Unsampled)
  x$stream_1_props <- dplyr::select(x$stream_1_props, !Unsampled)
  x
}
make_unsampled_signature <- function(delta = 0) {
  u_scaled <- colMeans(S_scaled) + delta * pc$rotation[, 1]
  u_original <- u_scaled * scale_tracer + center_tracer
  tibble::as_tibble_row(as.list(u_original)) |>
    dplyr::mutate(source = "Unsampled", .before = 1)
}
# adapting MixMustR function for this
run_mixmustr_models_sim <- function(model_choices, data_streams_list,
                                    tracer_list, sigma_ln_rho, ...) {
  mu_tab <- tracer_list$mus
  sig_tab <- tracer_list$sigmas
  n_tab <- tracer_list$ns
  sig_tab_required <- any(model_choices$sample_tracer | 
    (!model_choices$fix_unsampled & !model_choices$sample_tracer)
  )
  if (sig_tab_required) {
    sig_tab_er <- "You need a valid data.frame of tracers signature SDs."
    if (is.null(sig_tab)) {
      stop(sig_tab_er)
    } else {
      MixMustR:::check_sd_tabs(sig_tab, mu_tab)
    }
  }
  models <- vector(mode = "list", length = nrow(model_choices))
  for (i in seq_len(nrow(model_choices))) {
    # build_stancode(model_choices$sample_tracer[i],
    #                model_choices$fix_unsampled[i],
    #                model_choices$hierarchical[i],
    #                model_choices$code_path[i])
    model_path <- model_choices$code_path[i]
    sdata <- mixmustr_wrangle_input(
      data_streams_list, tracer_list, model_path, sigma_ln_rho = sigma_ln_rho,
      sample_tracer = model_choices$sample_tracer[i],
      fix_unsampled = model_choices$fix_unsampled[i],
      hierarchical = model_choices$hierarchical[i]
    )
    mod_name_suffix_i <- file_path_sans_ext(basename(model_path))
    timing <- system.time({
      model <- run_mixmod(model_path, mod_name_suffix_i, data = sdata, ...)
    })
    models[[i]] <- list(timing = timing, model = model, data = sdata)
    names(models)[i] <- mod_name_suffix_i
  }
  models
}
rename_sources <- function(df, source_names) {
  source_cols <- grep("^Source_[0-9]+$", names(df), value = TRUE)
  source_cols <- source_cols[
    order(as.integer(sub("^Source_", "", source_cols)))
  ]
  if (length(source_names) != length(source_cols)) {
    stop(
      "Length of `source_names` must match the number of Source_N columns.\n",
      "Got ", length(source_names), " source names but found ",
      length(source_cols), " Source_N columns: ",
      paste(source_cols, collapse = ", ")
    )
  }
  rename_map <- stats::setNames(source_cols, source_names)
  df |>
    dplyr::rename(dplyr::all_of(rename_map))
}
sim_fct <- function(prop_seed, tracer_seed, stan_seed, tracer_list, sign_df,
                    sim_id, model_id, truth_stream, to_write = TRUE,
                    source_names, noises, ...) {
  # 1. True proportions / stream agreement
  prop_dfs <- purrr::map(noises, function(noise_i, seed_i, ...) {
    do.call(
      simulate_mvn_mixture, c(list(...), list(seed = seed_i, noise = noise_i))
    ) |>
    (`[[`)("dataset") |>
    dplyr::mutate(group = paste0("G", group)) |>
    rename_sources(source_names)
  }, seed_i = prop_seed, ...)
  df_props_1 <- prop_dfs[[1]]
  df_props_2 <- prop_dfs[[2]]
  # 2. Tracer stochasticity
  bcs_si <- dplyr::select(sign_df, source:Study)
  if ("24:0" %in% names(sign_df)) {
    bcs_fa <- dplyr::select(sign_df, source, `24:0`:`20:5w3 (n)`, Study) |>
      dplyr::mutate(Taxa = NA) |>
      dplyr::relocate(Taxa, .after = source)
    synth_df <- MixMustR:::make_mixture_data(
      bcs_si, bcs_fa, df_props_1[, -1], df_props_2[, -1],
      truth_stream = truth_stream, synthetic_data, rand_gen = TRUE, sd_ = 1,
      seed = tracer_seed
    ) |>
      rm_unsampled()
  } else {
    synth_df <- make_mixture_data_si(
      bcs_si, df_props_1[, -1], df_props_2[, -1],
      truth_stream = truth_stream, synthetic_data, rand_gen = TRUE, sd_ = 1,
      seed = tracer_seed
    ) |>
      rm_unsampled()
  }
  # 3. MCMC reproducibility
  model_run <- run_mixmustr_models_sim(
    MixMustR::mixmustr_models[model_id, ], synth_df, tracer_list, sign_df,
    sigma_ln_rho = 1, iter = 2e3, warmup = 1e3, chains = 4, cores = 4,
    seed = stan_seed, control = list(adapt_delta = 0.99, max_treedepth = 15)
  )
  fit <- model_run[[1]]$model
  fit_standata <- model_run[[1]]$data
  if (to_write) {
    out_name <- paste0("revision_sim_", sim_id, ".RData")
    save(fit, fit_standata, synth_df, tracer_list, sign_df,
         file = file.path("copilot", "sim_models", out_name))
  }
  # Save and return diagnostics
  out <- diagnose_mixmustr_model_run(
    model_run = model_run[[1]],
    model_choice = MixMustR::mixmustr_models[model_id, ],
    synth_df = synth_df,
    tracer_list = tracer_list,
    sign_df = sign_df,
    sigma_ln_rho = 1,
    sim_id = sim_id,
    model_id = model_id,
    prop_seed = prop_seed,
    tracer_seed = tracer_seed,
    stan_seed = stan_seed,
    max_treedepth = 15L,
    output_dir = "copilot/diagnostics_rds"
  )
  # Drop heavy objects before return
  rm(fit, model_run)
  gc(verbose = FALSE)
  out
}
