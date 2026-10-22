reshape_tracer_pars <- function(tp, signatures) {
  # Unsampled means
  u_mu <- cbind(
    signatures$bcs_si[signatures$bcs_si$source == "Unsampled",
                      c("source", "d(13C/12C)", "d(15N/14N)")],
    signatures$bcs_fa[signatures$bcs_fa$source == "Unsampled",
                      c("24:0", "18:1w9", "18:2w6", "18:3w3",
                        "20:4w6", "20:5w3")]
  )
  # Unsampled SDs
  u_sigma <- data.frame(
    source = "Unsampled",
    `d(13C/12C)` = signatures$bcs_si[
      signatures$bcs_si$source == "Unsampled", "d(13C/12C) (SD)"
    ],
    `d(15N/14N)` = signatures$bcs_si[
      signatures$bcs_si$source == "Unsampled", "d(15N/14N) (SD)"
    ],
    `24:0` = signatures$bcs_fa[
      signatures$bcs_fa$source == "Unsampled", "24:0 (SD)"
    ],
    `18:1w9` = signatures$bcs_fa[
      signatures$bcs_fa$source == "Unsampled", "18:1w9 (SD)"
    ],
    `18:2w6` = signatures$bcs_fa[
      signatures$bcs_fa$source == "Unsampled", "18:2w6 (SD)"
    ],
    `18:3w3` = signatures$bcs_fa[
      signatures$bcs_fa$source == "Unsampled", "18:3w3 (SD)"
    ],
    `20:4w6` = signatures$bcs_fa[
      signatures$bcs_fa$source == "Unsampled", "20:4w6 (SD)"
    ],
    `20:5w3` = signatures$bcs_fa[
      signatures$bcs_fa$source == "Unsampled", "20:5w3 (SD)"
    ],
    check.names = FALSE
  )
  # n = 10 for every Unsampled tracer
  u_n <- u_sigma
  u_n[, -1] <- 10L
  # Append
  tp$mus    <- dplyr::bind_rows(tp$mus, u_mu)
  tp$sigmas <- dplyr::bind_rows(tp$sigmas, u_sigma)
  tp$ns     <- dplyr::bind_rows(tp$ns, u_n)
  # Check
  tp
}
regenarate_synth_df_iso <- function(tracer_parameters) {
  # Adapted from MixMustR source code
  mus <- tracer_parameters$mus
  source_names <- c(mus$source, "Unsampled")
  j <- length(source_names)
  synthetic_data <- data.frame(group = paste0("G", rep(1:10, each = 10)))
  signatures <- MixMustR:::augment_bcs_with_unsampled(
    signatures <- MixMustR::bcs_si, signatures <- MixMustR::bcs_fa, mus,
    delta = 3
  )
  # Stream 1 (tracers) and stream 2 (e.g. eDNA) share group centres (same seed)
  # but differ in within-group noise: low noise (0.5) vs low noise (1).
  df_props_1 <- MixMustR:::simulate_mvn_mixture(
    G = 10, M = 10, J = j, noise = 0.5, seed = 10004
  )$dataset |>
    dplyr::mutate(group = paste0("G", group)) |>
    MixMustR:::rename_sources(source_names)
  df_props_2 <- MixMustR:::simulate_mvn_mixture(
    G = 10, M = 10, J = j, noise = 1, seed = 10004
  )$dataset |>
    dplyr::mutate(group = paste0("G", group)) |>
    MixMustR:::rename_sources(source_names)
  list(
    df = MixMustR:::make_mixture_data(
      signatures$bcs_si, signatures$bcs_fa, df_props_1[, -1], df_props_2[, -1],
      truth_stream = 1, synthetic_data, rand_gen = TRUE, sd_ = 1, seed = 20004
    ),
    pars = reshape_tracer_pars(tracer_parameters, signatures)
  )
}
###########################
###########################
###########################
# MIXSIAR FUNCTIONS
###########################
###########################
###########################
generate_mixsiar_objts <- function(mu_tab, sig_tab, synth_df) {
  session_1 <- tempdir()
  sources_file <- file.path(session_1, "sources.csv")
  discr_file <- file.path(session_1, "discrimination.csv")
  mix_file <- file.path(session_1, "mixture.csv")
  # Write files first
  dplyr::left_join(
    mu_tab |>
      dplyr::select(!`24:0`:`20:5w3`) |>
      dplyr::rename(Meand13C = `d(13C/12C)`, Meand15N = `d(15N/14N)`),
    sig_tab |>
      dplyr::select(!`24:0`:`20:5w3`) |>
      dplyr::rename(SDd13C = `d(13C/12C)`, SDd15N = `d(15N/14N)`),
    by = dplyr::join_by(source)
  ) |>
    dplyr::arrange(source) |>
    dplyr::rename(X = source) |>
    dplyr::mutate(n = 10) |>
    write.csv(sources_file, row.names = FALSE)
  data.frame(
    X = sort(mu_tab$source), Meand13C = 0, SDd13C = 0, Meand15N = 0, SDd15N = 0
  ) |>
    write.csv(discr_file, row.names = FALSE)
  synth_df$df_stream_1 |>
    dplyr::select(random = group, d13C = `d(13C/12C)`, d15N = `d(15N/14N)`) |>
    write.csv(mix_file, row.names = FALSE)
  # Now load them back for MixSIAR purposes assuming group-level structure
  mix <- MixSIAR::load_mix_data(
    filename = mix_file, iso_names = c("d13C", "d15N"), factors = "random",
    fac_random = TRUE, fac_nested = NULL, cont_effects = NULL
  )
  source <- MixSIAR::load_source_data(
    filename = sources_file, source_factors = NULL, conc_dep = FALSE,
    data_type = "means", mix
  )
  discr <- MixSIAR::load_discr_data(filename = discr_file, mix)
  list(mix = mix, source = source, discr = discr)
}

plot_isospace <- function(isotopes, mix, source, discr, sim_df) {
  df <- data.frame(x = mix$data_iso[, isotopes[1]],
                   y = mix$data_iso[, isotopes[2]],
                   group = mix$data$random)
  x_label <- expression(paste(delta^13, "C (\u2030)", sep = ""))
  y_label <- expression(paste(delta^15, "N (\u2030)", sep = ""))
  index <- seq_len(source$n.sources)
  discr_mu_plot <- discr$mu
  discr_sig2_plot <- discr$sig2
  MU_plot <- array(NA, dim = c(length(source$S_MU[, 1]), 2))
  SIG_plot <- array(NA, dim = c(length(source$S_SIG[, 1]), 2))
  for (iso in 1:2) {
    MU_plot[, iso] <- source$S_MU[, isotopes[iso]] +
      discr_mu_plot[, isotopes[iso]]
    SIG_plot[, iso] <- sqrt(source$S_SIG[, isotopes[iso]]^2 +
      discr_sig2_plot[, isotopes[iso]])
  }
  df_sources <- data.frame(
    x = MU_plot[, 1], y = MU_plot[, 2], ymin = MU_plot[, 2] - SIG_plot[, 2],
    ymax = MU_plot[, 2] + SIG_plot[, 2], xmin = MU_plot[, 1] - SIG_plot[, 1],
    xmax = MU_plot[, 1] + SIG_plot[, 1],
    label = source$source_names
  )
  source_labels <- data.frame(
    x = MU_plot[index, 1] + c(-2, 0.5, -4.3, -4.1, -3, 0.5, -4.6),
    y = MU_plot[index, 2] + c(-0.4, -0.3, 0.3, -0.3, 0.3, 0.3, 0.3),
    label = gsub("Terrestrial", "T.", source$source_names)
  )
  sources <- sim_df |>
    dplyr::select(Macroalgae:`Unsampled`)
  out <- vector(mode = "list", length = nrow(sim_df))
  for (i in seq_along(out)) {
    if (i == 1) {
      y_pos <- 1
    } else {
      delta <- ifelse(sim_df$group[i] != sim_df$group[i - 1], 0.3, 0.1)
      y_pos <- y_pos - delta
    }
    x_start <- 0
    out_j <- vector(mode = "list", length = sum(sources[i, ] > 0))
    n_ <- 1
    for (j in seq_len(ncol(sources))) {
      if (sources[i, j] == 0) {
        next
      }
      xs <- c(x_start, x_start + sources[i, j])
      x_start <- xs[2]
      out_j[[n_]] <- data.frame(
        row = i, group = sim_df$group[i], source = names(sources)[j],
        y = y_pos, x = xs
      )
      n_ <- n_ + 1
    }
    out[[i]] <- dplyr::bind_rows(out_j)
  }
  out <- dplyr::bind_rows(out)
  ref_points <- out |>
    dplyr::group_by(group) |>
    dplyr::summarise(y = mean(y), x = -0.1) |>
    dplyr::ungroup()
  palcol <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
              "black", "#999999")
  shapes <- 1:10
  names(shapes) <- ref_points$group
  plot_1 <- ggplot(data = out) +
    geom_line(mapping = aes(x = x, y = y, group = row, colour = source),
              show.legend = FALSE) +
    geom_point(data = ref_points, mapping = aes(x = x, y = y, shape = group),
               show.legend = FALSE, size = 2.5) +
    annotate(geom = "text", x = -0.1, y = 1.5, hjust = 0.5, label = "Group") +
    scale_y_continuous(expand = c(0, 0), limits = c(-11.1, 1.7)) +
    scale_x_continuous(expand = c(0, 0), limits = c(-0.2, 1.01),
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_shape_manual(values = shapes) +
    scale_colour_manual(values = palcol) +
    labs(x = "Mixing proportion", y = NULL,
         title = "Simulated mixing proportions") +
    theme_bw() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 13),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank())
  plot_2 <- ggplot(data = df, mapping = aes(x = x, y = y)) +
    geom_point(mapping = aes(shape = group),
               show.legend = FALSE, size = 2.5) +
    scale_shape_manual(values = shapes) +
    geom_errorbar(
      data = df_sources, mapping = aes(ymin = ymin, ymax = ymax),
      linetype = 1, linewidth = 0.2, show.legend = FALSE, width = 0
    ) +
    geom_errorbar(
      data = df_sources, mapping = aes(xmin = xmin, xmax = xmax),
      linewidth = 0.2, width = 0, linetype = 1, show.legend = FALSE,
      orientation = "y"
    ) +
    geom_text(
      data = source_labels, mapping = aes(x = x, y = y, label = label),
      show.legend = FALSE, hjust = 0
    ) +
    geom_point(
      data = df_sources, mapping = aes(fill = label),
      size = 3, shape = 21, show.legend = FALSE
    ) +
    scale_fill_manual(values = palcol) +
    scale_x_continuous(expand = c(0, 0), limits = c(-29.5, -9)) +
    labs(x = x_label, y = y_label, title = "Isospace") +
    theme_bw() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 13),
          legend.position = "bottom",
          legend.text = element_text(size = 15),
          legend.title = element_blank(),
          legend.key.size = unit(1.5, "line")) +
    guides(shape = guide_legend(override.aes = list(size = 4)))
  plot_1 + plot_2 +
    patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")") &
    theme(plot.tag = element_text(face = "bold"))
}

prior_uncertainty_plot <- function(x) {
  ln_sigma_x <- rep(0.1, length(x))
  plot_a <- MixMustR::evaluate_uncertainty(x, ln_sigma_x)$plot +
    labs(
      title = substitute(italic(pi) == "{" * a * "}", list(a = paste0(x, collapse = ", "))),
      subtitle = substitute(italic(sigma)["log(" * italic("q") * ")"] == "{" * a
       *
       "}", list(a = paste0(ln_sigma_x, collapse = ", ")))
    ) +
    theme(plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 12))
  ln_sigma_x <- rep(1, length(x))
  plot_b <- MixMustR::evaluate_uncertainty(x, ln_sigma_x)$plot +
    labs(
      title = substitute(italic(pi) == "{" * a * "}", list(a = paste0(x, collapse = ", "))),
      subtitle = substitute(italic(sigma)["ln" * italic(rho)] == "{" * a * "}", list(a = paste0(ln_sigma_x, collapse = ", ")))
    ) +
    theme(plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 12))
  (plot_a + plot_b) +
    patchwork::plot_layout(axes = "collect") &
    theme(axis.title = element_text(size = 14))
}

###########################
###########################
###########################
# HPC SIMULATION FUNCTIONS
###########################
###########################
###########################

# -----------------------------------------------------------------------------
# User-facing entry point
# -----------------------------------------------------------------------------

summarise_hpc_simulations <- function(
  jobs_file             = "hpc/jobs_dev.csv",
  rds_dir               = "hpc/diagnostics_rds",
  rds_pattern           = "^diag_model_.*\\.rds$",
  out_dir               = "hpc/diagnostics_rds/summaries",
  rhat_threshold        = 1.01,
  bfmi_threshold        = 0.3,
  rmse_success_threshold = 0.05,
  sigma_ratio_lower     = 0.5,
  sigma_ratio_upper     = 2
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "tables_raw"),     recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "tables_summary"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "figures"),        recursive = TRUE, showWarnings = FALSE)

  jobs <- read_jobs_metadata(jobs_file)

  diag_files <- list.files(rds_dir, pattern = rds_pattern,
                           full.names = TRUE, recursive = FALSE)

  if (!length(diag_files)) {
    stop(
      "No diagnostic RDS files found in `rds_dir` ('", rds_dir, "') ",
      "matching pattern '", rds_pattern, "'.",
      call. = FALSE
    )
  }

  message("Reading ", length(diag_files), " diagnostic RDS file(s)...")
  diagnostics <- read_hpc_diagnostics(diag_files)

  message("Joining diagnostics to task-array metadata...")
  diagnostics <- join_all_diagnostics_to_jobs(diagnostics, jobs)

  message("Writing harmonised raw diagnostic tables...")
  write_diagnostic_tables(diagnostics, file.path(out_dir, "tables_raw"))

  message("Building manuscript-level summaries...")
  summaries <- build_manuscript_summaries(
    diagnostics            = diagnostics,
    rhat_threshold         = rhat_threshold,
    bfmi_threshold         = bfmi_threshold,
    rmse_success_threshold = rmse_success_threshold,
    sigma_ratio_lower      = sigma_ratio_lower,
    sigma_ratio_upper      = sigma_ratio_upper
  )

  write_diagnostic_tables(summaries, file.path(out_dir, "tables_summary"))

  message("Creating figures...")
  figs <- make_manuscript_figures(summaries, diagnostics)
  save_summary_figures(figs, file.path(out_dir, "figures"))

  message("Done. Outputs written to: ", out_dir)

  list(
    jobs        = jobs,
    diagnostics = diagnostics,
    summaries   = summaries,
    figures     = figs,
    files       = list(
      raw_tables     = list.files(file.path(out_dir, "tables_raw"),     full.names = TRUE),
      summary_tables = list.files(file.path(out_dir, "tables_summary"), full.names = TRUE),
      figures        = list.files(file.path(out_dir, "figures"),        full.names = TRUE)
    )
  )
}

# -----------------------------------------------------------------------------
# Metadata handling
# -----------------------------------------------------------------------------

read_jobs_metadata <- function(jobs_file) {
  jobs <- readr::read_csv(jobs_file, show_col_types = FALSE)

  # Fix 7: model_id added to the required-columns check.
  required <- c(
    "model_id", "tracers", "source_names", "degree", "sim_id", "G", "M",
    "input_file", "prop_seed", "tracer_seed", "stan_seed", "noises"
  )
  missing_cols <- setdiff(required, names(jobs))
  if (length(missing_cols)) {
    stop(
      "jobs file is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # truth_stream is expected in the hpc simulations; this fallback
  # keeps the summariser backwards-compatible with older task arrays.
  if (!"truth_stream" %in% names(jobs)) {
    warning(
      "`truth_stream` was not found in the jobs metadata; setting it to NA.",
      call. = FALSE
    )
    jobs$truth_stream <- NA_integer_
  }

  jobs |>
    mutate(
      tracers      = as.character(.data$tracers),
      source_names = as.character(.data$source_names),
      degree       = as.character(.data$degree),
      noises       = as.character(.data$noises),
      truth_stream = suppressWarnings(as.integer(.data$truth_stream)),
      truth_stream_label = dplyr::case_when(
        .data$truth_stream == 1L ~ "truth from stream 1",
        .data$truth_stream == 2L ~ "truth from stream 2",
        is.na(.data$truth_stream) ~ "truth stream unspecified",
        TRUE ~ paste0("truth from stream ", .data$truth_stream)
      ),
      n_tracers = count_semicolon_items(.data$tracers),
      n_sources = count_semicolon_items(.data$source_names),
      tracer_set = dplyr::case_when(
        .data$n_tracers == 2 ~ "2 tracers",
        TRUE ~ paste0(.data$n_tracers, " tracers")
      ),
      source_set = dplyr::case_when(
        .data$n_sources == 4 ~ "3 sampled + unsampled",
        .data$n_sources == 7 ~ "6 sampled + unsampled",
        TRUE ~ paste0(.data$n_sources - 1L, " sampled + unsampled")
      ),
      # Fix 6: dplyr::recode() replaced with dplyr::recode_values().
      stream_scenario = dplyr::recode_values(
        .data$degree,
        "con" ~ "convergent streams",
        "div" ~ "divergent streams",
        default = .data$degree
      ),
      noise_stream_1 = parse_noise(.data$noises, index = 1),
      noise_stream_2 = parse_noise(.data$noises, index = 2),
      condition_id = paste(
        .data$tracer_set,
        .data$source_set,
        .data$stream_scenario,
        .data$truth_stream_label,
        sep = " | "
      )
    )
}

count_semicolon_items <- function(x) {
  lengths(strsplit(as.character(x), ";", fixed = TRUE))
}

parse_noise <- function(x, index = 1) {
  map_dbl(strsplit(as.character(x), ";", fixed = TRUE), function(z) {
    z <- suppressWarnings(as.numeric(z))
    if (length(z) >= index) z[[index]] else NA_real_
  })
}

# -----------------------------------------------------------------------------
# RDS reading and harmonisation
# -----------------------------------------------------------------------------

read_hpc_diagnostics <- function(diag_files) {
  diag_list <- map(diag_files, read_one_hpc_rds)

  names_to_bind <- c(
    "stan_diagnostics", "prop_recovery", "true_resid_summary",
    "total_sigma", "mu_unsampled", "mu_source"
  )

  set_names(names_to_bind) |>
    map(function(nm) {
      map_dfr(diag_list, function(x) x[[nm]] %||% tibble())
    }) |>
    map(normalise_diagnostic_types)
}

read_one_hpc_rds <- function(path) {
  obj  <- readRDS(path)
  flat <- flatten_diagnostic_object(obj)

  if (!length(flat)) {
    warning("No diagnostic object recognised in ", basename(path), call. = FALSE)
    return(empty_diag_list(file = path))
  }

  out <- list(
    stan_diagnostics   = map_dfr(flat, pluck_or_empty, "stan_diagnostics"),
    prop_recovery      = map_dfr(flat, pluck_or_empty, "prop_recovery"),
    true_resid_summary = map_dfr(flat, pluck_or_empty, "true_resid_summary"),
    total_sigma        = map_dfr(flat, pluck_or_empty, "total_sigma"),
    mu_unsampled       = map_dfr(flat, pluck_or_empty, "mu_unsampled"),
    mu_source          = map_dfr(flat, pluck_or_empty, "mu_source")
  )

  map(out, ~ mutate(.x, diagnostic_file = basename(path), .before = 1))
}

empty_diag_list <- function(file = NA_character_) {
  nms <- c(
    "stan_diagnostics", "prop_recovery", "true_resid_summary",
    "total_sigma", "mu_unsampled", "mu_source"
  )
  set_names(nms) |>
    map(~ tibble(diagnostic_file = basename(file)))
}

flatten_diagnostic_object <- function(obj) {
  is_diag <- function(x) {
    is.list(x) && any(c(
      "stan_diagnostics", "prop_recovery", "true_resid_summary",
      "total_sigma", "mu_unsampled", "mu_source"
    ) %in% names(x))
  }

  if (is_diag(obj)) return(list(obj))

  if (is.list(obj)) {
    out <- list()
    recurse <- function(x) {
      if (is_diag(x)) {
        out[[length(out) + 1L]] <<- x
      } else if (is.list(x)) {
        walk(x, recurse)
      }
    }
    recurse(obj)
    return(out)
  }

  list()
}

pluck_or_empty <- function(x, nm) {
  y <- x[[nm]]
  if (is.null(y)) tibble() else as_tibble(y)
}

normalise_diagnostic_types <- function(x) {
  if (!nrow(x)) return(x)
  x |>
    mutate(
      across(any_of(c("sim_id", "model_id", "prop_seed", "tracer_seed",
                      "stan_seed", "truth_stream")), as.integer),
      across(any_of(c("sample_tracer", "fix_unsampled",
                      "hierarchical", "covered")), as.logical),
      across(any_of(c("variant", "target", "source", "tracer", "degree",
                      "truth_stream_label", "condition_id")), as.character)
    )
}

join_all_diagnostics_to_jobs <- function(diagnostics, jobs) {
  map(diagnostics, join_diagnostic_to_jobs, jobs = jobs)
}

join_diagnostic_to_jobs <- function(tbl, jobs) {
  if (!nrow(tbl)) return(tbl)

  join_cols <- intersect(
    c("sim_id", "prop_seed", "tracer_seed", "stan_seed"), names(tbl)
  )

  # sim_id is unique across all scenarios in jobs_dev.csv; if present, it
  # alone is sufficient and avoids accidental fan-out joins.
  by_cols <- if ("sim_id" %in% join_cols) "sim_id" else join_cols

  if (!length(by_cols)) {
    warning(
      "Could not find join columns in diagnostic table; ",
      "returning unjoined table.",
      call. = FALSE
    )
    return(tbl)
  }

  jobs_small <- jobs |>
    select(any_of(c(
      by_cols,
      "tracers", "source_names", "degree", "G", "M", "input_file", "noises",
      "n_tracers", "n_sources", "tracer_set", "source_set", "stream_scenario",
      "truth_stream", "truth_stream_label",
      "noise_stream_1", "noise_stream_2", "condition_id"
    ))) |>
    distinct()

  tbl |>
    left_join(jobs_small, by = by_cols, suffix = c("", "_job"))
}

write_diagnostic_tables <- function(x, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  iwalk(x, function(tbl, nm) {
    if (is.data.frame(tbl)) {
      readr::write_csv(tbl, file.path(out_dir, paste0(nm, ".csv")))
    }
  })
}

# -----------------------------------------------------------------------------
# Plotting helper
# -----------------------------------------------------------------------------

# Fix 4 + 8: safe facet_grid wrapper that:
#   * returns NULL (no facet) when all supplied variables are absent/empty;
#   * backtick-quotes each variable name so names containing parentheses or
#     other special formula characters do not silently corrupt the formula.
facet_grid_safe <- function(rows = NULL, cols = NULL, scales = "fixed") {
  rows <- rows[!is.na(rows) & nzchar(rows)]
  cols <- cols[!is.na(cols) & nzchar(cols)]
  if (!length(rows) && !length(cols)) return(NULL)

  btick    <- function(x) paste0("`", x, "`")
  row_expr <- if (length(rows)) paste(btick(rows), collapse = " + ") else "."
  col_expr <- if (length(cols)) paste(btick(cols), collapse = " + ") else "."

  ggplot2::facet_grid(
    stats::as.formula(paste(row_expr, "~", col_expr)),
    scales = scales
  )
}

# -----------------------------------------------------------------------------
# Manuscript-level summaries
# -----------------------------------------------------------------------------

build_manuscript_summaries <- function(
  diagnostics,
  rhat_threshold         = 1.01,
  bfmi_threshold         = 0.3,
  rmse_success_threshold = 0.05,
  sigma_ratio_lower      = 0.5,
  sigma_ratio_upper      = 2
) {
  list(
    simulation_design = summarise_simulation_design(
      diagnostics$stan_diagnostics
    ),
    stan_by_model_condition = summarise_stan_by_model_condition(
      diagnostics$stan_diagnostics,
      rhat_threshold = rhat_threshold,
      bfmi_threshold = bfmi_threshold
    ),
    prop_recovery_by_model_condition = summarise_prop_recovery_by_model_condition(
      diagnostics$prop_recovery,
      rmse_success_threshold = rmse_success_threshold
    ),
    prop_recovery_by_source = summarise_prop_recovery_by_source(
      diagnostics$prop_recovery
    ),
    sigma_calibration_by_model_condition = summarise_sigma_calibration_by_model_condition(
      diagnostics$total_sigma,
      sigma_ratio_lower = sigma_ratio_lower,
      sigma_ratio_upper = sigma_ratio_upper
    ),
    sigma_calibration_by_tracer = summarise_sigma_calibration_by_tracer(
      diagnostics$total_sigma
    ),
    mu_unsampled_recovery = summarise_signature_recovery(
      diagnostics$mu_unsampled,
      grouping_cols = c(
        "model_id", "variant", "condition_id", "tracer_set", "source_set",
        "stream_scenario", "truth_stream", "truth_stream_label", "tracer"
      ),
      prefix = "mu_unsampled"
    ),
    mu_source_recovery = summarise_signature_recovery(
      diagnostics$mu_source,
      grouping_cols = c(
        "model_id", "variant", "condition_id", "tracer_set", "source_set",
        "stream_scenario", "truth_stream", "truth_stream_label", "source", "tracer"
      ),
      prefix = "mu_source"
    ),
    model_ranking_by_condition = rank_models_by_condition(
      prop_recovery    = diagnostics$prop_recovery,
      stan_diagnostics = diagnostics$stan_diagnostics,
      rhat_threshold   = rhat_threshold
    ),
    truth_stream_contrast = summarise_truth_stream_contrast(
      diagnostics$prop_recovery
    )
  )
}

summarise_simulation_design <- function(stan_tbl) {
  if (!nrow(stan_tbl)) return(tibble())

  stan_tbl |>
    distinct(
      sim_id, model_id, variant, sample_tracer, fix_unsampled, hierarchical,
      condition_id, tracer_set, source_set, stream_scenario,
      truth_stream, truth_stream_label,
      n_tracers, n_sources, degree, G, M, noises,
      prop_seed, tracer_seed, stan_seed
    ) |>
    summarise(
      n_model_assessments     = n(),
      n_unique_simulations    = n_distinct(sim_id),
      n_models                = n_distinct(model_id),
      n_conditions            = n_distinct(condition_id),
      n_tracer_sets           = n_distinct(tracer_set),
      n_source_sets           = n_distinct(source_set),
      n_stream_scenarios      = n_distinct(stream_scenario),
      n_truth_streams         = n_distinct(truth_stream, na.rm = TRUE),
      .groups = "drop"
    )
}

summarise_stan_by_model_condition <- function(stan_tbl,
                                              rhat_threshold = 1.01,
                                              bfmi_threshold = 0.3) {
  if (!nrow(stan_tbl)) return(tibble())

  # Fix 2: guard optional column that may be absent in older diagnostic runs.
  has_params_excluded <-
    "n_params_excluded_from_diagnostics" %in% names(stan_tbl)

  stan_tbl |>
    group_by(
      model_id, variant, sample_tracer, fix_unsampled, hierarchical,
      condition_id, tracer_set, source_set, stream_scenario,
      truth_stream, truth_stream_label
    ) |>
    summarise(
      n_fits                       = n(),
      total_divergences            = sum(n_divergent,      na.rm = TRUE),
      prop_fits_with_divergences   = mean(n_divergent > 0, na.rm = TRUE),
      total_max_treedepth          = sum(n_max_treedepth,      na.rm = TRUE),
      prop_fits_with_treedepth     = mean(n_max_treedepth > 0, na.rm = TRUE),
      median_max_rhat              = median(max_rhat,           na.rm = TRUE),
      q95_max_rhat                 = quantile(max_rhat, 0.95,   na.rm = TRUE),
      prop_fits_rhat_gt_threshold  = mean(max_rhat > rhat_threshold, na.rm = TRUE),
      # Fix 5: replace Inf (returned when all values are NA) with NA_real_.
      min_bfmi = {
        v <- min(min_bfmi, na.rm = TRUE)
        if (is.finite(v)) v else NA_real_
      },
      prop_fits_bfmi_below_threshold = mean(min_bfmi < bfmi_threshold, na.rm = TRUE),
      median_min_n_eff               = median(min_n_eff, na.rm = TRUE),
      median_params_excluded = if (has_params_excluded) {
        median(n_params_excluded_from_diagnostics, na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) |>
    arrange(condition_id, model_id)
}

summarise_prop_recovery_by_model_condition <- function(prop_tbl,
                                                       rmse_success_threshold = 0.05) {
  if (!nrow(prop_tbl)) return(tibble())
  prop_tbl |>
    dplyr::summarise(
      n_source_summaries = dplyr::n(),
      median_rmse = median(rmse, na.rm = TRUE),
      mean_rmse = mean(rmse, na.rm = TRUE),
      q05_rmse = quantile(rmse, 0.05, na.rm = TRUE),
      q95_rmse = quantile(rmse, 0.95, na.rm = TRUE),
      median_mae = median(mae, na.rm = TRUE),
      mean_mae = mean(mae, na.rm = TRUE),
      q05_mae = quantile(mae, 0.05, na.rm = TRUE),
      q95_mae = quantile(mae, 0.95, na.rm = TRUE),
      median_abs_bias = median(abs(bias), na.rm = TRUE),
      mean_abs_bias = mean(abs(bias), na.rm = TRUE),
      q05_abs_bias = quantile(abs(bias), 0.05, na.rm = TRUE),
      q95_abs_bias = quantile(abs(bias), 0.95, na.rm = TRUE),
      median_cor = median(cor, na.rm = TRUE),
      mean_cor = mean(cor, na.rm = TRUE),
      q05_cor = quantile(cor, 0.05, na.rm = TRUE),
      q95_cor = quantile(cor, 0.95, na.rm = TRUE),
      min_cor = min(cor, na.rm = TRUE),
      max_cor = max(cor, na.rm = TRUE),
      prop_source_summaries_rmse_below_threshold = mean(
        rmse < rmse_success_threshold, na.rm = TRUE
      ),
      .by = c(
        model_id, variant, sample_tracer, fix_unsampled, hierarchical,
        condition_id, tracer_set, source_set, stream_scenario, truth_stream,
        truth_stream_label, target
      )
    ) |>
    dplyr::arrange(condition_id, target, median_rmse)
}

summarise_prop_recovery_by_source <- function(prop_tbl) {
  if (!nrow(prop_tbl)) return(tibble())
  prop_tbl |>
    dplyr::summarise(
      n = dplyr::n(),
      median_rmse = median(rmse, na.rm = TRUE),
      mean_rmse = mean(rmse, na.rm = TRUE),
      q05_rmse = quantile(rmse, 0.05, na.rm = TRUE),
      q95_rmse = quantile(rmse, 0.95, na.rm = TRUE),
      median_mae = median(mae, na.rm = TRUE),
      mean_mae = mean(mae, na.rm = TRUE),
      q05_mae = quantile(mae, 0.05, na.rm = TRUE),
      q95_mae = quantile(mae, 0.95, na.rm = TRUE),
      median_abs_bias = median(abs(bias), na.rm = TRUE),
      mean_abs_bias = mean(abs(bias), na.rm = TRUE),
      q05_abs_bias = quantile(abs(bias), 0.05, na.rm = TRUE),
      q95_abs_bias = quantile(abs(bias), 0.95, na.rm = TRUE),
      median_cor = median(cor, na.rm = TRUE),
      mean_cor = mean(cor, na.rm = TRUE),
      q05_cor = quantile(cor, 0.05, na.rm = TRUE),
      q95_cor = quantile(cor, 0.95, na.rm = TRUE),
      min_cor = min(cor, na.rm = TRUE),
      max_cor = max(cor, na.rm = TRUE),
      .by = c(
        model_id, variant, condition_id, tracer_set, source_set,
        stream_scenario, truth_stream, truth_stream_label, target, source
      )
    ) |>
    dplyr::arrange(condition_id, target, model_id, desc(median_rmse))
}

summarise_sigma_calibration_by_model_condition <- function(total_sigma_tbl,
                                                           sigma_ratio_lower = 0.5,
                                                           sigma_ratio_upper = 2) {
  if (!nrow(total_sigma_tbl)) return(tibble())
  if (!all(c("ratio", "abs_log_ratio") %in% names(total_sigma_tbl))) return(tibble())

  total_sigma_tbl |>
    group_by(
      model_id, variant, sample_tracer, fix_unsampled, hierarchical,
      condition_id, tracer_set, source_set, stream_scenario,
      truth_stream, truth_stream_label
    ) |>
    summarise(
      n_tracer_summaries       = n(),
      median_sigma_ratio       = median(ratio,          na.rm = TRUE),
      mean_sigma_ratio         = mean(ratio,            na.rm = TRUE),
      q05_sigma_ratio          = quantile(ratio, 0.05,  na.rm = TRUE),
      q95_sigma_ratio          = quantile(ratio, 0.95,  na.rm = TRUE),
      median_abs_log_ratio     = median(abs_log_ratio,  na.rm = TRUE),
      mean_abs_log_ratio       = mean(abs_log_ratio,    na.rm = TRUE),
      prop_ratios_within_bounds = mean(
        ratio >= sigma_ratio_lower & ratio <= sigma_ratio_upper, na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    arrange(condition_id, model_id)
}

summarise_sigma_calibration_by_tracer <- function(total_sigma_tbl) {
  if (!nrow(total_sigma_tbl)) return(tibble())
  if (!all(c("ratio", "abs_log_ratio") %in% names(total_sigma_tbl))) return(tibble())

  total_sigma_tbl |>
    group_by(
      model_id, variant, condition_id, tracer_set, source_set,
      stream_scenario, truth_stream, truth_stream_label, tracer
    ) |>
    summarise(
      n                    = n(),
      median_sigma_ratio   = median(ratio,          na.rm = TRUE),
      mean_sigma_ratio     = mean(ratio,            na.rm = TRUE),
      q05_sigma_ratio      = quantile(ratio, 0.05,  na.rm = TRUE),
      q95_sigma_ratio      = quantile(ratio, 0.95,  na.rm = TRUE),
      median_abs_log_ratio = median(abs_log_ratio,  na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(condition_id, model_id, desc(median_abs_log_ratio))
}

summarise_signature_recovery <- function(tbl, grouping_cols,
                                         prefix = "signature") {
  if (!nrow(tbl)) return(tibble())
  if (!all(c("covered", "abs_error_median") %in% names(tbl))) return(tibble())

  tbl |>
    group_by(across(all_of(grouping_cols))) |>
    summarise(
      n                  = n(),
      coverage           = mean(covered,            na.rm = TRUE),
      median_abs_error   = median(abs_error_median, na.rm = TRUE),
      mean_abs_error     = mean(abs_error_median,   na.rm = TRUE),
      q95_abs_error      = quantile(abs_error_median, 0.95, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(condition_id, desc(median_abs_error)) |>
    mutate(signature_type = prefix, .before = 1)
}

rank_models_by_condition <- function(prop_recovery, stan_diagnostics,
                                     rhat_threshold = 1.01) {
  if (!nrow(prop_recovery) || !nrow(stan_diagnostics)) return(tibble())

  prop_rank <- prop_recovery |>
    filter(target == "stream_1_props") |>
    group_by(
      model_id, variant, sample_tracer, fix_unsampled, hierarchical,
      condition_id, tracer_set, source_set, stream_scenario,
      truth_stream, truth_stream_label
    ) |>
    summarise(
      median_rmse    = median(rmse,      na.rm = TRUE),
      mean_rmse      = mean(rmse,        na.rm = TRUE),
      median_cor     = median(cor,       na.rm = TRUE),
      mean_abs_bias  = mean(abs(bias),   na.rm = TRUE),
      .groups = "drop"
    )

  stan_rank <- stan_diagnostics |>
    group_by(model_id, variant, condition_id, truth_stream, truth_stream_label) |>
    summarise(
      total_divergences   = sum(n_divergent,        na.rm = TRUE),
      total_max_treedepth = sum(n_max_treedepth,    na.rm = TRUE),
      prop_rhat_flag      = mean(max_rhat > rhat_threshold, na.rm = TRUE),
      .groups = "drop"
    )

  prop_rank |>
    left_join(
      stan_rank,
      by = c("model_id", "variant", "condition_id",
             "truth_stream", "truth_stream_label")
    ) |>
    mutate(
      diagnostic_penalty =
        10L * as.integer(total_divergences > 0) +
        2L  * as.integer(total_max_treedepth > 0) +
        5   * prop_rhat_flag,
      manuscript_score = median_rmse + 0.01 * diagnostic_penalty
    ) |>
    group_by(condition_id) |>
    arrange(manuscript_score, desc(median_cor), .by_group = TRUE) |>
    mutate(rank_within_condition = row_number()) |>
    ungroup()
}

summarise_truth_stream_contrast <- function(prop_tbl) {
  if (!nrow(prop_tbl)) return(tibble())
  if (!"truth_stream" %in% names(prop_tbl)) return(tibble())

  prop_tbl |>
    filter(target == "stream_1_props") |>
    group_by(
      model_id, variant, sample_tracer, fix_unsampled, hierarchical,
      tracer_set, source_set, stream_scenario, truth_stream, truth_stream_label
    ) |>
    summarise(
      median_rmse   = median(rmse,      na.rm = TRUE),
      mean_rmse     = mean(rmse,        na.rm = TRUE),
      median_cor    = median(cor,       na.rm = TRUE),
      mean_abs_bias = mean(abs(bias),   na.rm = TRUE),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      id_cols     = c(model_id, variant, sample_tracer, fix_unsampled,
                      hierarchical, tracer_set, source_set, stream_scenario),
      names_from  = truth_stream_label,
      values_from = c(median_rmse, mean_rmse, median_cor, mean_abs_bias),
      names_sep   = "__"
    )
}

# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------

make_manuscript_figures <- function(summaries, diagnostics) {
  list(
    prop_rmse_by_condition       = plot_prop_rmse_by_condition(summaries$prop_recovery_by_model_condition),
    stan_flags_by_condition      = plot_stan_flags(summaries$stan_by_model_condition),
    sigma_calibration_by_condition = plot_sigma_calibration(summaries$sigma_calibration_by_model_condition),
    unsampled_signature_coverage = plot_signature_coverage(summaries$mu_unsampled_recovery),
    model_ranking                = plot_model_ranking(summaries$model_ranking_by_condition)
  )
}

# Fix 1 + 4: truth_stream_label is now retained through the pivot and
# all facets use facet_grid_safe() so missing/all-NA columns don't crash.
plot_stan_flags <- function(tbl) {
  if (!nrow(tbl)) return(NULL)

  long <- tbl |>
    select(
      model_id, condition_id, tracer_set, source_set,
      stream_scenario, truth_stream_label,
      prop_fits_with_divergences, prop_fits_with_treedepth,
      prop_fits_rhat_gt_threshold
    ) |>
    tidyr::pivot_longer(
      tidyselect::starts_with("prop_fits"),
      names_to  = "diagnostic",
      values_to = "proportion"
    )

  ggplot(data = long,
         mapping = aes(x = as.factor(model_id), y = proportion,
                       fill = diagnostic)) +
    geom_col(position = position_dodge(width = 0.8)) +
    facet_grid_safe(rows = c("tracer_set", "source_set"),
                    cols = c("stream_scenario", "truth_stream_label")) +
    labs(
      x     = "Model variant",
      y     = "Proportion of fits flagged",
      fill  = "Diagnostic flag",
      title = "Stan diagnostic flags across simulation conditions"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_prop_rmse_by_condition <- function(tbl) {
  if (!nrow(tbl)) return(NULL)

  ggplot(
    data = dplyr::filter(tbl, target == "stream_1_props"),
    mapping = aes(
      x    = forcats::fct_reorder(as.factor(model_id), median_rmse),
      y    = median_rmse,
      fill = stream_scenario
    )
  ) +
    geom_col(position = position_dodge(width = 0.8)) +
    facet_grid_safe(rows = c("tracer_set", "truth_stream_label"),
                    cols = c("source_set"), scales = "free_x") +
    labs(
      x     = "Model variant",
      y     = "Median RMSE for stream-1 proportions",
      fill  = "Stream scenario",
      title = "Mixing-proportion recovery across model variants and simulation conditions"
    ) +
    theme_bw()
}

plot_sigma_calibration <- function(tbl) {
  if (!nrow(tbl)) return(NULL)

  ggplot(data = tbl,
         mapping = aes(x = as.factor(model_id), y = median_sigma_ratio,
                       colour = stream_scenario)) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    facet_grid_safe(rows = c("tracer_set", "truth_stream_label"),
                    cols = c("source_set")) +
    labs(
      x      = "Model variant",
      y      = "Median fitted/empirical residual SD ratio",
      colour = "Stream scenario",
      title  = "Residual-scale calibration across model variants"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_signature_coverage <- function(tbl) {
  if (!nrow(tbl)) return(NULL)

  ggplot(data = tbl,
         mapping = aes(x = tracer, y = coverage,
                       colour = as.factor(model_id))) +
    geom_hline(yintercept = 0.95, linetype = 2) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    facet_grid_safe(rows = c("tracer_set", "source_set"),
                    cols = c("stream_scenario", "truth_stream_label"),
                    scales = "free_x") +
    labs(
      x      = "Tracer",
      y      = "Coverage of true unsampled signature",
      colour = "Model variant",
      title  = "Unsampled-source signature recovery"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_model_ranking <- function(tbl) {
  if (!nrow(tbl)) return(NULL)

  ggplot(
    data = dplyr::filter(tbl, rank_within_condition <= 8),
    mapping = aes(
      x    = reorder(as.factor(model_id), manuscript_score),
      y    = manuscript_score,
      fill = stream_scenario
    )
  ) +
    geom_col(position = position_dodge(width = 0.8)) +
    facet_grid_safe(rows = c("tracer_set", "truth_stream_label"),
                    cols = c("source_set"), scales = "free_x") +
    labs(
      x     = "Model variant",
      y     = "Composite manuscript score",
      fill  = "Stream scenario",
      title = "Model ranking by proportion recovery and Stan diagnostics"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

save_summary_figures <- function(figs, out_dir, width = 12, height = 8,
                                 dpi = 300) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  iwalk(figs, function(p, nm) {
    if (!is.null(p)) {
      ggplot2::ggsave(
        filename = file.path(out_dir, paste0(nm, ".png")),
        plot     = p,
        width    = width,
        height   = height,
        dpi      = dpi
      )
    }
  })
}

# -----------------------------------------------------------------------------
# Optional manuscript helper: compact text-ready tables
# -----------------------------------------------------------------------------

make_text_ready_summary <- function(summaries) {
  list(
    best_models_by_condition = summaries$model_ranking_by_condition |>
      filter(rank_within_condition == 1) |>
      select(
        condition_id, tracer_set, source_set, stream_scenario,
        truth_stream, truth_stream_label,
        model_id, sample_tracer, fix_unsampled, hierarchical,
        median_rmse, median_cor, mean_abs_bias,
        total_divergences, total_max_treedepth, prop_rhat_flag,
        manuscript_score
      ),

    overall_model_performance = summaries$model_ranking_by_condition |>
      group_by(model_id, variant, sample_tracer, fix_unsampled, hierarchical) |>
      summarise(
        median_rmse_across_conditions  = median(median_rmse,          na.rm = TRUE),
        mean_rmse_across_conditions    = mean(mean_rmse,              na.rm = TRUE),
        median_cor_across_conditions   = median(median_cor,           na.rm = TRUE),
        mean_abs_bias_across_conditions = mean(mean_abs_bias,         na.rm = TRUE),
        median_rank                    = median(rank_within_condition, na.rm = TRUE),
        n_best_conditions              = sum(rank_within_condition == 1, na.rm = TRUE),
        .groups = "drop"
      ) |>
      arrange(median_rank, median_rmse_across_conditions)
  )
}

# -----------------------------------------------------------------------------
# Theme and panel labels
# -----------------------------------------------------------------------------

theme_redesign <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(
        fill = "grey90", colour = "grey35", linewidth = 0.35
      ),
      strip.text = element_text(face = "bold", size = base_size + 0.75),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title = element_text(face = "bold"),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(size = base_size),
      plot.margin = margin(2, 2, 2, 2)
    )
}

add_panel_label <- function(plot, label, title, subtitle = NULL) {
  title_block <- cowplot::ggdraw() +
    cowplot::draw_label(
      paste0(label, ")"),
      x = 0.01, y = 0.68, hjust = 0, vjust = 0.5, fontface = "bold", size = 18
    ) +
    cowplot::draw_label(
      title, x = 0.075, y = 0.82, hjust = 0, vjust = 0.5, fontface = "bold",
      size = 13.5
    )
  if (!is.null(subtitle)) {
    title_block <- title_block +
      cowplot::draw_label(
        subtitle, x = 0.075, y = 0.44, hjust = 0, vjust = 0.5, size = 10.8
      )
  }
  title_block / plot + patchwork::plot_layout(heights = c(0.075, 0.925))
}

# -----------------------------------------------------------------------------
# Data labelling helpers
# -----------------------------------------------------------------------------
make_model_variant_labels <- function(dat) {
  dat |>
    mutate(
      tracer_part = dplyr::if_else(sample_tracer, "S", "F"),
      unsampled_part = dplyr::if_else(fix_unsampled, "F", "S"),
      hierarchy_part = dplyr::if_else(hierarchical, "H", "I"),
      model_key = paste0(tracer_part, unsampled_part, "-", hierarchy_part),
      model_key = factor(
        model_key,
        levels = c(
          "FF-H", "FF-I", "FS-H", "FS-I", "SF-H", "SF-I", "SS-H", "SS-I"
        )
      ),
      model_long = paste0(
        model_key, " = ",
        dplyr::if_else(sample_tracer, "sampled tracer", "fixed tracer"), " + ",
        dplyr::if_else(fix_unsampled, "fixed unsampled", "sampled unsampled"),
        " + ", dplyr::if_else(hierarchical, "hierarchical", "independent")
      )
    )
}

clean_simulation_labels <- function(dat) {
  want <- c("stream_scenario", "truth_stream_label", "tracer_set", "source_set")
  if (!all(want %in% names(dat))) {
    stop("Wanted columns not present.")
  }
  dat |>
    dplyr::mutate(
      degree_label = dplyr::case_when(
        stringr::str_detect(
          stream_scenario, stringr::regex("convergent", ignore_case = TRUE)
        ) ~ "Convergent",
        stringr::str_detect(
          stream_scenario, stringr::regex("divergent", ignore_case = TRUE)
        ) ~ "Divergent",
        .default = as.character(stream_scenario)
      ),
      # Normalise truth_stream_label to the capitalised "Truth stream N" form
      # expected by the colour scales in plot_detail_with_separator(),
      # regardless of the casing produced by the upstream summariser (e.g.
      # "truth from stream 1").
      truth_stream_label = dplyr::case_when(
        grepl("1", as.character(truth_stream_label), fixed = TRUE) ~ "Truth stream 1",
        grepl("2", as.character(truth_stream_label), fixed = TRUE) ~ "Truth stream 2",
        .default = as.character(truth_stream_label)
      ),
      tracer_set_label = as.character(tracer_set),
      source_set_label = as.character(source_set),
      condition_short = paste(degree_label, tracer_set_label, sep = "\n") |>
        factor(
          levels = c(
            "Convergent\n2 tracers", "Convergent\n8 tracers",
            "Divergent\n2 tracers", "Divergent\n8 tracers"
          )
        ),
      source = forcats::fct_relevel(
        source,
        "Macroalgae", "Mangrove", "Saltmarsh", "Seagrass", "SPOM",
        "Terrestrial grass", "Unsampled"
      )
    )
}

# -----------------------------------------------------------------------------
# Compact horizontal colour keys
# -----------------------------------------------------------------------------

make_horizontal_key <- function(title, low, high, mid = NULL, midpoint = 0,
                                limits = c(0, 1), breaks = waiver(),
                                labels = waiver(), bar_height = 0.02) {
  key_dat <- tibble::tibble(
    x = seq(limits[1], limits[2], length.out = 256),
    y = 1
  )
  p <- ggplot(data = key_dat, mapping = aes(x = x, y = y, fill = x)) +
    geom_tile(height = bar_height) +
    annotate(
      geom = "rect", xmin = limits[1], xmax = limits[2],
      ymin = 1 - bar_height / 2, ymax = 1 + bar_height / 2,
      fill = NA, colour = "black", linewidth = 0.22
    ) +
    scale_x_continuous(breaks = breaks, labels = labels, expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = title) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(
        face = "bold", size = 9, hjust = 0.5,
        margin = margin(b = 5)
      ),
      axis.text.x = element_text(size = 8.5, margin = margin(t = 1)),
      plot.margin = margin(0, 2, 0, 2)
    )
  if (is.null(mid)) {
    p + scale_fill_gradient(
      low = low, high = high, limits = limits,
      oob = scales::squish, guide = "none"
    )
  } else {
    p + scale_fill_gradient2(
      low = low, mid = mid, high = high,
      midpoint = midpoint, limits = limits,
      oob = scales::squish, guide = "none"
    )
  }
}

wrap_key <- function(key_plot) {
  # Keeps each colour bar narrow and vertically centred in its row.
  ((patchwork::plot_spacer() | key_plot | patchwork::plot_spacer()) +
     patchwork::plot_layout(widths = c(0.18, 0.64, 0.18))) /
    patchwork::plot_spacer() +
    patchwork::plot_layout(heights = c(0.36, 0.64))
}

make_model_key_plot <- function() {
  # Compact two-column legend with v12-like horizontal gap to avoid overlap.
  cowplot::ggdraw() +
    cowplot::draw_label(
      "Model variant key",
      x = 0.55, y = 0.93,
      hjust = 0.5, vjust = 1,
      fontface = "bold", size = 11.8
    ) +
    cowplot::draw_label(
      "FF = fixed tracer + fixed unsampled",
      x = 0, y = 0.62,
      hjust = 0, vjust = 0.5,
      size = 10.3
    ) +
    cowplot::draw_label(
      "FS = fixed tracer + sampled unsampled",
      x = 0, y = 0.39,
      hjust = 0, vjust = 0.5,
      size = 10.3
    ) +
    cowplot::draw_label(
      "SF = sampled tracer + fixed unsampled",
      x = 0, y = 0.16,
      hjust = 0, vjust = 0.5,
      size = 10.3
    ) +
    cowplot::draw_label(
      "SS = sampled tracer + sampled unsampled",
      x = 0.65, y = 0.62,
      hjust = 0, vjust = 0.5,
      size = 10.3
    ) +
    cowplot::draw_label(
      "H  = hierarchical",
      x = 0.65, y = 0.39,
      hjust = 0, vjust = 0.5,
      size = 10.3
    ) +
    cowplot::draw_label(
      "I  = independent",
      x = 0.65, y = 0.16,
      hjust = 0, vjust = 0.5,
      size = 10.3
    )
}

make_panel_a_key_column <- function(key_rmse, key_mae, key_bias, key_cor) {
  # Bars are compact and centred on the y-midpoint of each heatmap row.
  # Positions are calibrated for heatmap spacer rows of 0.17.
  cowplot::ggdraw() +
    cowplot::draw_plot(
      key_rmse, x = 0.18, y = 0.833, width = 0.64, height = 0.112
    ) +
    cowplot::draw_plot(
      key_mae,  x = 0.18, y = 0.575, width = 0.64, height = 0.112
    ) +
    cowplot::draw_plot(
      key_bias, x = 0.18, y = 0.317, width = 0.64, height = 0.112
    ) +
    cowplot::draw_plot(
      key_cor,  x = 0.18, y = 0.059, width = 0.64, height = 0.112
    )
}

# -----------------------------------------------------------------------------
# Panel A
# -----------------------------------------------------------------------------

prepare_panel_a_data <- function(df, target_to_plot = "stream_1_props") {
  df |>
    dplyr::filter(target == .env$target_to_plot) |>
    make_model_variant_labels() |>
    clean_simulation_labels() |>
    dplyr::select(
      sim_id, model_id, model_key, degree_label, tracer_set_label,
      condition_short, source_set_label, truth_stream_label, source, rmse, mae,
      bias, cor
    ) |>
    tidyr::pivot_longer(
      cols = c("rmse", "mae", "bias", "cor"), names_to = "metric",
      values_to = "value"
    ) |>
    dplyr::mutate(
      metric = factor(
        metric, levels = c("rmse", "mae", "bias", "cor"),
        labels = c("RMSE", "MAE", "Bias", "Correlation")
      )
    ) |>
    # Global summary: collapse truth stream and source-set richness.
    dplyr::summarise(
      value = median(value, na.rm = TRUE),
      .by = c(metric, condition_short, model_key, source)
    )
}

make_metric_heatmap <- function(dat, metric_label, scale_fill,
                                show_x_strips = FALSE, show_x_axis = FALSE) {
  p <- dat |>
    dplyr::filter(metric == .env$metric_label) |>
    ggplot(mapping = aes(x = source, y = model_key, fill = value)) +
      geom_tile(colour = "white", linewidth = 0.25) +
      facet_grid(. ~ condition_short, scales = "free_x", space = "free_x") +
      scale_fill +
      labs(x = NULL, y = metric_label) +
      theme_redesign(base_size = 10) +
      theme(
        axis.text.y = element_text(size = 10),
        axis.text.x = if (show_x_axis) {
          element_text(size = 9, angle = 45, hjust = 1)
        } else {
          element_blank()
        },
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.title.y = element_text(
          size = 12, face = "bold", angle = 90,
          vjust = 0.5, margin = margin(r = 5)
        ),
        legend.position = "none",
        panel.spacing.x = grid::unit(0.12, "cm"),
        plot.margin = margin(1, 1, 1, 1)
      )
  if (!show_x_strips) {
    p <- p +
      theme(
        strip.text.x = element_blank(),
        strip.background.x = element_blank()
      )
  }
  p
}

make_panel_a <- function(panel_a_dat, rmse_limits = c(0, 0.4),
                         mae_limits = c(0, 0.3), bias_limits = c(-0.25, 0.25),
                         cor_limits = c(0, 1)) {
  p_a_rmse <- make_metric_heatmap(
    panel_a_dat, "RMSE",
    scale_fill_gradient(
      low = "#fff7bc", high = "#d7191c",
      limits = rmse_limits, oob = scales::squish
    ),
    show_x_strips = TRUE, show_x_axis = FALSE
  )
  p_a_mae <- make_metric_heatmap(
    panel_a_dat, "MAE",
    scale_fill_gradient(
      low = "#fff7bc", high = "#d7191c",
      limits = mae_limits, oob = scales::squish
    ),
    show_x_strips = FALSE, show_x_axis = FALSE
  )
  p_a_bias <- make_metric_heatmap(
    panel_a_dat, "Bias",
    scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#b2182b",
      midpoint = 0, limits = bias_limits, oob = scales::squish
    ),
    show_x_strips = FALSE, show_x_axis = FALSE
  )
  p_a_cor <- make_metric_heatmap(
    panel_a_dat, "Correlation",
    scale_fill_gradient(
      low = "#f7fcb9", high = "#238443",
      limits = cor_limits, oob = scales::squish
    ),
    show_x_strips = FALSE, show_x_axis = TRUE
  )
  heatmap_stack <- p_a_rmse /
    patchwork::plot_spacer() /
    p_a_mae /
    patchwork::plot_spacer() /
    p_a_bias /
    patchwork::plot_spacer() /
    p_a_cor +
    patchwork::plot_layout(heights = c(1, -0.1, 1, -0.1, 1, -0.1, 1))
  key_rmse <- make_horizontal_key(
    title = "RMSE", low = "#fff7bc", high = "#d7191c",
    limits = rmse_limits, breaks = pretty(rmse_limits, n = 3),
    labels = pretty(rmse_limits, n = 3), bar_height = 0.02
  )
  key_mae <- make_horizontal_key(
    title = "MAE", low = "#fff7bc", high = "#d7191c",
    limits = mae_limits, breaks = pretty(mae_limits, n = 3),
    labels = pretty(mae_limits, n = 3), bar_height = 0.02
  )
  key_bias <- make_horizontal_key(
    title = "Bias", low = "#2166ac", mid = "white", high = "#b2182b",
    midpoint = 0, limits = bias_limits,
    breaks = c(bias_limits[1], 0, bias_limits[2]),
    labels = c(as.character(bias_limits[1]), "0", as.character(bias_limits[2])),
    bar_height = 0.02
  )
  key_cor <- make_horizontal_key(
    title = "Correlation", low = "#f7fcb9", high = "#238443",
    limits = cor_limits, breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1"), bar_height = 0.02
  )
  key_col <- make_panel_a_key_column(
    key_rmse = key_rmse, key_mae = key_mae, key_bias = key_bias,
    key_cor = key_cor
  )
  top_row <- cowplot::plot_grid(
    heatmap_stack, key_col, nrow = 1, rel_widths = c(5.25, 1.0), align = "h",
    axis = "tb"
  )
  # Put the model-variant key below the whole vertical stack and centred across
  # all x axes, with an explicit spacer row of about 0.5 cm above it.
  model_key_row <- cowplot::ggdraw() +
    # Heatmap region occupies about 84% of Panel A width. This x/width centres
    # the compact two-column legend under the middle of all heatmap x axes.
    cowplot::draw_plot(
      make_model_key_plot(), x = 0.155, y = 0.00, width = 0.52, height = 1
    )
  panel_a_inner <- cowplot::plot_grid(
    top_row, NULL, model_key_row, ncol = 1,
    # Reduce the gap above the key and allocate more height to the heatmaps.
    rel_heights = c(0.815, 0.020, 0.165)
  )
  add_panel_label(
    panel_a_inner, label = "A", title = "Global performance summary",
    subtitle = "Performance across convergence states, tracer sets, model variants, and sources"
  )
}

# -----------------------------------------------------------------------------
# Panel B
# -----------------------------------------------------------------------------

prepare_detail_axis <- function(dat) {
  out <- dat |>
    dplyr::mutate(source = as.character(source)) |>
    dplyr::group_by(detail_block) |>
    dplyr::mutate(
      source_id = as.integer(factor(source, levels = unique(source)))
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      block_id = as.integer(detail_block),
      x = source_id +
        if_else(block_id == 2L, max(source_id, na.rm = TRUE) + 1.2, 0),
      x_lab = source
    )
  separator_x <- max(out$x[out$block_id == 1], na.rm = TRUE) + 0.6
  block_labs <- out |>
    dplyr::summarise(x_mid = mean(range(x)), .by = c(detail_block, block_id))
  x_labs <- out |>
    dplyr::distinct(x, x_lab) |>
    dplyr::arrange(x)
  list(data = out, separator_x = separator_x, block_labs = block_labs,
       x_labs = x_labs)
}

plot_detail_with_separator <- function(dat, title) {
  truth_cols <- c("Truth stream 1" = "#ef6c5b", "Truth stream 2" = "#39b7c4")
  shape_values <- c("2 tracers" = 21, "8 tracers" = 22)
  fill_values <- c("3 sampled sources" = "white",
                   "6 sampled sources" = "grey25")
  prep <- prepare_detail_axis(dat)
  dat2 <- prep$data |>
    dplyr::mutate(
      dodge_group = interaction(
        truth_stream_label, tracer_shape, source_fullness, drop = TRUE
      ),
      source_fullness = factor(
        source_fullness, levels = c("3 sampled sources", "6 sampled sources")
      )
    )
  pd <- position_dodge(width = 0.45)
  y_top <- max(dat2$upper, na.rm = TRUE)
  y_label <- y_top * 1.2
  ggplot(
    data = dat2,
    mapping = aes(x = x, y = mean_rmse, colour = truth_stream_label,
                  shape = tracer_shape, fill = source_fullness,
                  group = dodge_group)
  ) +
    geom_vline(
      xintercept = prep$separator_x, linetype = 2, colour = "grey55",
      linewidth = 0.35
    ) +
    geom_errorbar(
      mapping = aes(ymin = lower, ymax = upper),
      width = 0, linewidth = 0.35, position = pd
    ) +
    geom_point(size = 2.8, stroke = 0.85, position = pd) +
    geom_text(
      data = prep$block_labs,
      mapping = aes(x = x_mid, y = y_label, label = detail_block),
      inherit.aes = FALSE, size = 4.3, lineheight = 0.9, vjust = 0.9
    ) +
    scale_x_continuous(
      breaks = prep$x_labs$x, labels = prep$x_labs$x_lab,
      expand = expansion(mult = c(0.04, 0.04))
    ) +
    scale_y_continuous(
      limits = c(0, y_label * 1.06), expand = expansion(mult = c(0, 0.02))
    ) +
    scale_colour_manual(values = truth_cols, name = "Truth stream") +
    scale_shape_manual(values = shape_values, name = "Tracer set") +
    scale_fill_manual(values = fill_values, name = "Sampled-source set") +
    labs(x = NULL, y = "RMSE", title = title) +
    theme_redesign(base_size = 10) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1, size = 11),
      axis.text.y = element_text(size = 11),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.margin = margin(2, 3, 2, 3)
    )
}

make_panel_b_legend <- function() {
  truth_cols <- c("Truth stream 1" = "#ef6c5b", "Truth stream 2" = "#39b7c4")
  # Three rows with vertical separation. The legend remains custom-drawn.
  # So Panel B keeps the exact same layout and stays self-contained.
  y_truth <- 0.86
  y_tracer <- 0.49
  y_source <- 0.12
  cowplot::ggdraw() +
    cowplot::draw_label(
      "Truth stream:", x = 0.02, y = y_truth, hjust = 0, fontface = "bold",
      size = 10.5
    ) +
    cowplot::draw_grob(
      grid::pointsGrob(
        x = grid::unit(0.36, "npc"), y = grid::unit(y_truth, "npc"), pch = 16,
        gp = grid::gpar(col = truth_cols[["Truth stream 1"]], fill = truth_cols[["Truth stream 1"]], cex = 0.95, lwd = 1)
      )
    ) +
    cowplot::draw_label("1", x = 0.4, y = y_truth, hjust = 0, size = 10) +
    cowplot::draw_grob(
      grid::pointsGrob(
        x = grid::unit(0.49, "npc"), y = grid::unit(y_truth, "npc"), pch = 16,
        gp = grid::gpar(col = truth_cols[["Truth stream 2"]], fill = truth_cols[["Truth stream 2"]], cex = 0.95, lwd = 1)
      )
    ) +
    cowplot::draw_label("2", x = 0.53, y = y_truth, hjust = 0, size = 10) +
    cowplot::draw_label(
      "N. tracers:", x = 0.02, y = y_tracer, hjust = 0, fontface = "bold",
      size = 10.5
    ) +
    cowplot::draw_grob(
      grid::pointsGrob(
        x = grid::unit(0.3, "npc"), y = grid::unit(y_tracer, "npc"), pch = 21,
        gp = grid::gpar(col = "grey20", fill = "white", cex = 0.95, lwd = 1)
      )
    ) +
    cowplot::draw_label("2", x = 0.34, y = y_tracer, hjust = 0, size = 10) +
    cowplot::draw_grob(
      grid::pointsGrob(
        x = grid::unit(0.43, "npc"), y = grid::unit(y_tracer, "npc"), pch = 22,
        gp = grid::gpar(col = "grey20", fill = "white", cex = 0.95, lwd = 1)
      )
    ) +
    cowplot::draw_label("8", x = 0.47, y = y_tracer, hjust = 0, size = 10) +
    cowplot::draw_label(
      "N. sampled sources:", x = 0.02, y = y_source, hjust = 0,
      fontface = "bold", size = 10.5
    ) +
    cowplot::draw_grob(
      grid::pointsGrob(
        x = grid::unit(0.5, "npc"), y = grid::unit(y_source, "npc"), pch = 21,
        gp = grid::gpar(col = "grey20", fill = "white", cex = 0.95, lwd = 1)
      )
    ) +
    cowplot::draw_label("3", x = 0.54, y = y_source, hjust = 0, size = 10) +
    cowplot::draw_grob(
      grid::pointsGrob(
        x = grid::unit(0.63, "npc"), y = grid::unit(y_source, "npc"), pch = 21,
        gp = grid::gpar(col = "grey20", fill = "grey25", cex = 0.95, lwd = 1)
      )
    ) +
    cowplot::draw_label("6", x = 0.67, y = y_source, hjust = 0, size = 10)
}

make_panel_b <- function(prop_raw, target_to_plot = "stream_1_props") {
  detail_dat <- prop_raw |>
    dplyr::filter(target == .env$target_to_plot) |>
    make_model_variant_labels() |>
    clean_simulation_labels() |>
    dplyr::mutate(
      tracer_shape = dplyr::case_when(
        tracer_set_label == "2 tracers" ~ "2 tracers",
        tracer_set_label == "8 tracers" ~ "8 tracers",
        .default = tracer_set_label
      ),
      source_fullness = dplyr::case_when(
        source_set_label == "3 sampled + unsampled" ~ "3 sampled sources",
        source_set_label == "6 sampled + unsampled" ~ "6 sampled sources",
        .default = source_set_label
      )
    )
  best_3_models <- detail_dat |>
    dplyr::summarise(
      median_rmse = median(rmse, na.rm = TRUE),
      .by = c(model_id, model_key)
    ) |>
    dplyr::arrange(median_rmse) |>
    dplyr::slice_head(n = 3)
  best_sources <- c("Macroalgae", "Seagrass", "SPOM")
  best_detail_sum <- detail_dat |>
    dplyr::semi_join(best_3_models, by = c("model_id", "model_key")) |>
    dplyr::filter(
      source %in% .env$best_sources,
      (degree_label == "Convergent" & tracer_set_label == "2 tracers") |
        (degree_label == "Divergent" & tracer_set_label == "8 tracers")
    ) |>
    dplyr::mutate(
      detail_block = dplyr::case_when(
        degree_label == "Convergent" & tracer_set_label == "2 tracers" ~ "Convergent\n2 tracers",
        degree_label == "Divergent" & tracer_set_label == "8 tracers" ~ "Divergent\n8 tracers"
      ),
      detail_block = factor(
        detail_block,
        levels = c("Convergent\n2 tracers", "Divergent\n8 tracers")
      )
    ) |>
    dplyr::summarise(
      mean_rmse = mean(rmse, na.rm = TRUE),
      sd_rmse = sd(rmse, na.rm = TRUE),
      .by = c(detail_block, source, truth_stream_label, tracer_shape,
              source_fullness)
    ) |>
    dplyr::mutate(
      lower = pmax(mean_rmse - sd_rmse, 0), upper = mean_rmse + sd_rmse
    )
  hardest_sources <- detail_dat |>
    dplyr::summarise(
      median_rmse = median(rmse, na.rm = TRUE), .by = source
    ) |>
    dplyr::arrange(dplyr::desc(median_rmse)) |>
    dplyr::slice_head(n = 2) |>
    dplyr::pull(source) |>
    as.character()
  hard_detail_sum <- detail_dat |>
    dplyr::filter(
      source %in% .env$hardest_sources,
      (degree_label == "Convergent" & tracer_set_label == "8 tracers") |
        (degree_label == "Divergent" & tracer_set_label == "2 tracers")
    ) |>
    dplyr::mutate(
      detail_block = dplyr::case_when(
        degree_label == "Convergent" & tracer_set_label == "8 tracers" ~ "Convergent\n8 tracers",
        degree_label == "Divergent" & tracer_set_label == "2 tracers" ~ "Divergent\n2 tracers"
      ),
      detail_block = factor(
        detail_block,
        levels = c("Convergent\n8 tracers", "Divergent\n2 tracers")
      )
    ) |>
    dplyr::summarise(
      mean_rmse = mean(rmse, na.rm = TRUE),
      sd_rmse = sd(rmse, na.rm = TRUE),
      .by = c(detail_block, source, truth_stream_label, tracer_shape,
              source_fullness)
    ) |>
    dplyr::mutate(
      lower = pmax(mean_rmse - sd_rmse, 0), upper = mean_rmse + sd_rmse
    )
  p_b_best <- plot_detail_with_separator(
    best_detail_sum, "Best 3 model variants"
  )
  p_b_hard <- plot_detail_with_separator(hard_detail_sum, "Hardest sources")
  # Add an explicit spacer between the two detail plots to make the separation
  # more manuscript-friendly in the narrow Panel B column.
  panel_b_inner <- p_b_best / patchwork::plot_spacer() / p_b_hard / make_panel_b_legend() +
    patchwork::plot_layout(heights = c(1, 0.11, 1, 0.34))
  add_panel_label(
    patchwork::wrap_elements(panel_b_inner),
    label = "B",
    title = "   Focused detail (examples)",
    subtitle = "   Uncertainty and truth-stream contrasts"
  )
}

# -----------------------------------------------------------------------------
# User-facing functions
# -----------------------------------------------------------------------------

make_perform_figs <- function(results,
                                             target_to_plot = "stream_1_props",
                                             rmse_limits = c(0, 0.4),
                                             mae_limits = c(0, 0.3),
                                             bias_limits = c(-0.25, 0.25),
                                             cor_limits = c(0, 1),
                                             truth_stream = NULL) {
  prop_raw <- results$diagnostics$prop_recovery
  if (is.null(truth_stream)) {
    prop_raw_a <- prop_raw
  } else {
    prop_raw_a <- dplyr::filter(prop_raw, truth_stream == .env$truth_stream)
  }
  panel_a_dat <- prepare_panel_a_data(
    prop_raw_a, target_to_plot = target_to_plot
  )
  panel_a <- make_panel_a(
    panel_a_dat, rmse_limits = rmse_limits, mae_limits = mae_limits,
    bias_limits = bias_limits, cor_limits = cor_limits
  )
  panel_b <- make_panel_b(prop_raw, target_to_plot = target_to_plot)
  body <- cowplot::plot_grid(
    panel_a, panel_b, nrow = 1, rel_widths = c(2, 1), align = "h",
    axis = "tb"
  )
  cowplot::ggdraw() +
    cowplot::draw_label(
      "MixMustR simulation performance results",
      x = 0.5, y = 0.985, hjust = 0.5, vjust = 1, fontface = "bold", size = 20
    ) +
    # Leave a visible bottom margin for captions/figure spacing in manuscripts.
    cowplot::draw_plot(body, x = 0, y = 0.055, width = 1, height = 0.895)
}

###########################
###########################
###########################
# STAN DIAGNOSTIC PLOTS
###########################
###########################
###########################

prepare_stan_diagnostics_for_plotting <- function(results) {
  results$diagnostics$stan_diagnostics |>
    dplyr::mutate(
      model_id = factor(
        model_id,
        levels = 1:8,
        labels = c(
          "1 - SF-H", "2 - FF-H", "3 - SS-H", "4 - FS-H",
          "5 - SF-I", "6 - FF-I", "7 - SS-I", "8 - FS-I"
        )
      ),
      stream_scenario = factor(
        stream_scenario,
        levels = c("convergent streams", "divergent streams")
      ),
      truth_stream_label = factor(
        truth_stream_label,
        levels = c("truth from stream 1", "truth from stream 2")
      )
    )
}

# Remove range-skewing outliers separately for each diagnostic using
# Tukey's 1.5-IQR rule. For zero-inflated count diagnostics, where IQR = 0,
# retain observations up to the 99.5th percentile instead.
trim_diagnostic <- function(data, variable, discrete = FALSE) {
  x <- data[[variable]]
  q <- stats::quantile(x, c(0.25, 0.75), na.rm = TRUE, names = FALSE)
  iqr <- diff(q)

  if (iqr > 0) {
    lower <- q[1] - 1.5 * iqr
    upper <- q[2] + 1.5 * iqr
  } else {
    lower <- min(x, na.rm = TRUE)
    upper <- stats::quantile(x, 0.995, na.rm = TRUE, names = FALSE)
  }

  if (discrete) lower <- max(0, floor(lower))

  dplyr::filter(
    data,
    !is.na(.data[[variable]]),
    .data[[variable]] >= lower,
    .data[[variable]] <= upper
  )
}

make_ridgeline <- function(data, variable, figure_label, x_label,
                           trim = FALSE, ...) {
  plot_data <- trim_diagnostic(data, variable, discrete = FALSE)

  model_colours <- stats::setNames(
    viridisLite::viridis(nlevels(plot_data$model_id), option = "D"),
    levels(plot_data$model_id)
  )

  p <- ggplot(
    plot_data,
    aes(
      x = .data[[variable]],
      y = forcats::fct_rev(model_id),
      fill = model_id
    )
  )

  # Binned ridgelines are used for integer-valued diagnostics because these
  # distributions contain large point masses at zero. Continuous diagnostics
  # use kernel-density ridgelines.
  p <- p +
    ggridges::geom_density_ridges(
      scale = 1.15,
      rel_min_height = 0.001,
      alpha = 0.85,
      colour = "white",
      linewidth = 0.2
    )
  p <- p +
    facet_grid(
      rows = vars(truth_stream_label),
      cols = vars(stream_scenario),
      scales = "free_x"
    ) +
    scale_fill_manual(values = model_colours, guide = "none") +
    labs(
      title = figure_label,
      x = x_label,
      y = "Model variant"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      strip.background = element_rect(fill = "grey95", colour = "grey75"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 8),
      panel.spacing = grid::unit(0.6, "lines")
    )
  if (trim) {
    p + scale_x_continuous(...)
  } else {
    p
  }
}

plot_stan_diagnostic <- function(results, variable, figure_label, x_label,
                                 trim = FALSE, ...) {
  make_ridgeline(
    prepare_stan_diagnostics_for_plotting(results), variable, figure_label,
    x_label, trim = trim, ...
  )
}
