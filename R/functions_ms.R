###########################
###########################
###########################
# POST-PROCESSING FUNCTIONS
###########################
###########################
###########################
wrangle_r2 <- function(..., rm = TRUE) {
  x <- MixMustR::mixmustr_bayes_R2(...)
  if (rm) {
    x <- dplyr::select(x, !Source)
  }
  x
}

make_post_rand_long <- function(modfit, mu_tab, synth_df) {
  reff_df <- synth_df$sim_data |>
    dplyr::mutate(R = as.numeric(as.factor(paste0(Location, Habitat)))) |>
    dplyr::distinct(R, Location, Habitat) |>
    dplyr::arrange(R)
  dplyr::left_join(    
    (summary(modfit)$summary) |>
      data.frame() |>
      tibble::rownames_to_column("Variable") |>
      dplyr::filter(grepl("randvec\\[", Variable)) |>
      tidyr::separate_wider_regex(
        Variable, c("randvec\\[", R = "\\d+", "\\,", source = "\\d+", "\\]")
      ) |>
      dplyr::mutate(dplyr::across(c(R, source), as.numeric)) |>
      dplyr::select(mean, R, source) |>
      tidyr::pivot_wider(names_from = "source", values_from = "mean") |>
      MixMustR:::assign_new_names(c("R", mu_tab$Group, "Unsampled")) |>
      dplyr::arrange(R) |>
      dplyr::left_join(reff_df, by = "R") |>
      tidyr::pivot_longer(
        cols = !c(R, Location, Habitat), names_to = "source",
        values_to = "Predicted"
      ) |>
      dplyr::mutate(Predicted = exp(Predicted)),
    reshape_ref_data(synth_df, order_ref = mu_tab$Group) |>
      abs_log() |>
      data.frame(check.names = FALSE) |>
      (\(x, y)cbind(x, y))(
        x = dplyr::select(synth_df$props_df, Location, Habitat)
      ) |>
      dplyr::mutate(R = as.numeric(as.factor(paste0(Location, Habitat)))) |>
      tidyr::pivot_longer(
        cols = !c(R, Location, Habitat), names_to = "source",
        values_to = "Observed"
      ) |>
      dplyr::mutate(Observed = exp(Observed)) |>
      dplyr::group_by(R, Location, Habitat, source) |>
      dplyr::summarise(Observed = mean(Observed)) |>
      dplyr::ungroup(),
    by = join_by(R, Location, Habitat, source)
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
    x = MU_plot[index, 1] + c(-4, 0.5, 0.5, -4.3, -4.1, 0.5),
    y = MU_plot[index, 2] + c(-0.4, 0.3, 0.3, -0.3, -0.3, 0.3),
    label = gsub("Terrestrial", "T.", source$source_names)
  )
  sources <- sim_df |>
    dplyr::select(Macroalgae:`Terrestrial grass`)
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
  palcol <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#999999")
  shapes <- c(1:6, 8)
  names(shapes) <- ref_points$group
  plot_1 <- ggplot(data = out) +
    geom_line(mapping = aes(x = x, y = y, group = row, colour = source),
              show.legend = FALSE) +
    geom_point(data = ref_points, mapping = aes(x = x, y = y, shape = group),
               show.legend = FALSE, size = 2.5) +
    annotate(geom = "text", x = -0.1, y = 1.5, hjust = 0.5, label = "Group") +
    scale_y_continuous(expand = c(0, 0), limits = c(-8.1, 1.7)) +
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
    geom_errorbarh(
      data = df_sources, mapping = aes(xmin = xmin, xmax = xmax),
      linewidth = 0.2, height = 0, linetype = 1, show.legend = FALSE
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

# Modified from `MixSIAR::run_model`
mixsiar_run <- function (run, mix, source, discr, model_filename,
                        alpha.prior = 1, resid_err = NULL, process_err = NULL) {
  err_raw <- read.table(model_filename, comment.char = "", 
      sep = ":", skip = 7, nrows = 1, colClasses = "character")
  if (err_raw[1, 2] == " Residual only") 
      err <- "resid"
  if (err_raw[1, 2] == " Process only (MixSIR, for N = 1)") 
      err <- "process"
  if (err_raw[1, 2] == " Residual * Process") 
      err <- "mult"
  if (length(alpha.prior) == 1) {
      if (alpha.prior == 1) 
          alpha.prior = rep(1, source$n.sources)
  }
  if (!is.numeric(alpha.prior)) {
      stop(paste("*** Error: Your prior is not a numeric vector of length(n.sources).\n        Try again or choose the uninformative prior option. For example,\n        c(1,1,1,1) is a valid (uninformative) prior for 4 sources. ***", 
          sep = ""))
  }
  if (length(alpha.prior) != source$n.sources) {
      stop(paste("*** Error: Length of your prior does not match the\n        number of sources (", 
          source$n.sources, "). Try again. ***", sep = ""))
  }
  if (length(which(alpha.prior == 0)) != 0) {
      stop(paste("*** Error: You cannot set any alpha = 0.\n      Instead, set = 0.01.***", 
          sep = ""))
  }
  if (is.numeric(alpha.prior) == F) 
      alpha.prior = 1
  if (length(alpha.prior) == 1) 
      alpha = rep(alpha.prior, source$n.sources)
  if (length(alpha.prior) > 1 & length(alpha.prior) != source$n.sources) 
      alpha = rep(1, source$n.sources)
  if (length(alpha.prior) > 1 & length(alpha.prior) == source$n.sources) 
      alpha = alpha.prior
  if (is.list(run)) {
      mcmc <- run
  }
  else {
      if (run == "test") 
          mcmc <- list(chainLength = 1000, burn = 500, thin = 1, 
              chains = 3, calcDIC = TRUE)
      if (run == "very short") 
          mcmc <- list(chainLength = 10000, burn = 5000, thin = 5, 
              chains = 3, calcDIC = TRUE)
      if (run == "short") 
          mcmc <- list(chainLength = 50000, burn = 25000, thin = 25, 
              chains = 3, calcDIC = TRUE)
      if (run == "normal") 
          mcmc <- list(chainLength = 1e+05, burn = 50000, thin = 50, 
              chains = 3, calcDIC = TRUE)
      if (run == "long") 
          mcmc <- list(chainLength = 3e+05, burn = 2e+05, thin = 100, 
              chains = 3, calcDIC = TRUE)
      if (run == "very long") 
          mcmc <- list(chainLength = 1e+06, burn = 5e+05, thin = 500, 
              chains = 3, calcDIC = TRUE)
      if (run == "extreme") 
          mcmc <- list(chainLength = 3e+06, burn = 1500000, 
              thin = 500, chains = 3, calcDIC = TRUE)
  }
  n.sources <- source$n.sources
  N <- mix$N
  e <- matrix(rep(0, n.sources * (n.sources - 1)), nrow = n.sources, 
      ncol = (n.sources - 1))
  for (i in 1:(n.sources - 1)) {
      e[, i] <- exp(c(rep(sqrt(1/(i * (i + 1))), i), -sqrt(i/(i + 
          1)), rep(0, n.sources - i - 1)))
      e[, i] <- e[, i]/sum(e[, i])
  }
  cross <- array(data = NA, dim = c(N, n.sources, n.sources - 
      1))
  tmp.p <- array(data = NA, dim = c(N, n.sources))
  jags.params <- c("p.global", "loglik", "p.ind")
  f.data <- character(0)
  if (mix$n.effects > 0 & !mix$fere) {
      factor1_levels <- mix$FAC[[1]]$levels
      Factor.1 <- mix$FAC[[1]]$values
      cross.fac1 <- array(data = NA, dim = c(factor1_levels, 
          n.sources, n.sources - 1))
      tmp.p.fac1 <- array(data = NA, dim = c(factor1_levels, 
          n.sources))
      if (mix$FAC[[1]]$re) 
          jags.params <- c(jags.params, "p.fac1", "ilr.fac1", 
              "fac1.sig")
      else jags.params <- c(jags.params, "p.fac1", "ilr.fac1")
      f.data <- c(f.data, "factor1_levels", "Factor.1", "cross.fac1", 
          "tmp.p.fac1")
  }
  if (mix$n.effects > 1 & !mix$fere) {
      factor2_levels <- mix$FAC[[2]]$levels
      Factor.2 <- mix$FAC[[2]]$values
      if (mix$fac_nested[1]) {
          factor2_lookup <- mix$FAC[[1]]$lookup
          f.data <- c(f.data, "factor2_lookup")
      }
      if (mix$fac_nested[2]) {
          factor1_lookup <- mix$FAC[[2]]$lookup
          f.data <- c(f.data, "factor1_lookup")
      }
      cross.fac2 <- array(data = NA, dim = c(factor2_levels, 
          n.sources, n.sources - 1))
      tmp.p.fac2 <- array(data = NA, dim = c(factor2_levels, 
          n.sources))
      if (mix$FAC[[2]]$re) 
          jags.params <- c(jags.params, "p.fac2", "ilr.fac2", 
              "fac2.sig")
      else jags.params <- c(jags.params, "p.fac2", "ilr.fac2")
      f.data <- c(f.data, "factor2_levels", "Factor.2", "cross.fac2", 
          "tmp.p.fac2")
  }
  if (mix$fere) {
      factor1_levels <- mix$FAC[[1]]$levels
      Factor.1 <- mix$FAC[[1]]$values
      if (mix$n.re == 1) {
          cross.fac1 <- array(data = NA, dim = c(factor1_levels, 
              n.sources, n.sources - 1))
          tmp.p.fac1 <- array(data = NA, dim = c(factor1_levels, 
              n.sources))
          jags.params <- c(jags.params, "p.fac1")
          f.data <- c(f.data, "cross.fac1", "tmp.p.fac1")
      }
      factor2_levels <- mix$FAC[[2]]$levels
      Factor.2 <- mix$FAC[[2]]$values
      if (mix$FAC[[2]]$re) 
          jags.params <- c(jags.params, "fac2.sig")
      jags.params <- c(jags.params, "ilr.global", "ilr.fac1", 
          "ilr.fac2")
      f.data <- c(f.data, "factor1_levels", "Factor.1", "factor2_levels", 
          "Factor.2")
  }
  if (source$data_type == "raw") {
      SOURCE_array <- source$SOURCE_array
      n.rep <- source$n.rep
      s.data <- c("SOURCE_array", "n.rep")
  }
  else {
      MU_array <- source$MU_array
      SIG2_array <- source$SIG2_array
      n_array <- source$n_array
      s.data <- c("MU_array", "SIG2_array", "n_array")
  }
  if (!is.na(source$by_factor)) {
      source_factor_levels <- source$S_factor_levels
      s.data <- c(s.data, "source_factor_levels")
  }
  if (source$conc_dep) {
      conc <- source$conc
      s.data <- c(s.data, "conc")
  }
  c.data <- rep(NA, mix$n.ce)
  if (mix$n.ce > 0) {
      for (ce in 1:mix$n.ce) {
          name <- paste("Cont.", ce, sep = "")
          assign(name, as.vector(mix$CE[[ce]]))
          c.data[ce] <- paste("Cont.", ce, sep = "")
          jags.params <- c(jags.params, "ilr.global", paste("ilr.cont", 
              ce, sep = ""), "p.ind")
      }
  }
  X_iso <- mix$data_iso
  n.iso <- mix$n.iso
  frac_mu <- discr$mu
  frac_sig2 <- discr$sig2
  all.data <- c("X_iso", "N", "n.sources", "n.iso", "alpha", 
      "frac_mu", "e", "cross", "tmp.p")
  jags.data <- c(all.data, f.data, s.data, c.data)
  I <- diag(n.iso)
  if (err == "resid" && mix$n.iso > 1) 
      jags.data <- c(jags.data, "I")
  if (err != "resid") 
      jags.data <- c(jags.data, "frac_sig2")
  if (err == "mult") 
      jags.params <- c(jags.params, "resid.prop")
  jags.inits <- function() {
      list(p.global = as.vector(MCMCpack::rdirichlet(1, alpha)))
  }
  for (j in 1:n.iso) {
      if (source$data_type == "raw") {
          if (!is.na(source$by_factor)) {
              mean.pool <- mean(c(X_iso[, j], as.vector(SOURCE_array[, 
                j, , ])), na.rm = T)
              sd.pool <- sd(c(X_iso[, j], as.vector(SOURCE_array[, 
                j, , ])), na.rm = T)
              SOURCE_array[, j, , ] <- (SOURCE_array[, j, , 
                ] - mean.pool)/sd.pool
          }
          else {
              mean.pool <- mean(c(X_iso[, j], as.vector(SOURCE_array[, 
                j, ])), na.rm = T)
              sd.pool <- sd(c(X_iso[, j], as.vector(SOURCE_array[, 
                j, ])), na.rm = T)
              SOURCE_array[, j, ] <- (SOURCE_array[, j, ] - 
                mean.pool)/sd.pool
          }
      }
      else {
          if (!is.na(source$by_factor)) {
              mean.pool <- (N * mean(X_iso[, j], na.rm = T) + 
                as.vector(as.vector(n_array) %*% as.vector(MU_array[, 
                  j, ])))/sum(c(as.vector(n_array), N))
              if (N > 1) 
                sd.pool <- sqrt((sum((as.vector(n_array) - 
                  1) * as.vector(SIG2_array[, j, ])) + as.vector(as.vector(n_array) %*% 
                  as.vector(MU_array[, j, ])^2) + (N - 1) * 
                  stats::var(X_iso[, j], na.rm = T) + N * mean(X_iso[, 
                  j], na.rm = T)^2 - sum(c(as.vector(n_array), 
                  N)) * mean.pool^2)/(sum(c(as.vector(n_array), 
                  N)) - 1))
              if (N == 1) 
                sd.pool <- sqrt((sum((as.vector(n_array) - 
                  1) * as.vector(SIG2_array[, j, ])) + as.vector(as.vector(n_array) %*% 
                  as.vector(MU_array[, j, ])^2) + N * mean(X_iso[, 
                  j], na.rm = T)^2 - sum(c(as.vector(n_array), 
                  N)) * mean.pool^2)/(sum(c(as.vector(n_array), 
                  N)) - 1))
              MU_array[, j, ] <- (MU_array[, j, ] - mean.pool)/sd.pool
              SIG2_array[, j, ] <- SIG2_array[, j, ]/sd.pool^2
          }
          else {
              mean.pool <- (N * mean(X_iso[, j], na.rm = T) + 
                as.vector(as.vector(n_array) %*% as.vector(MU_array[, 
                  j])))/sum(c(as.vector(n_array), N))
              if (N > 1) 
                sd.pool <- sqrt((sum((as.vector(n_array) - 
                  1) * as.vector(SIG2_array[, j])) + as.vector(as.vector(n_array) %*% 
                  as.vector(MU_array[, j])^2) + (N - 1) * stats::var(X_iso[, 
                  j], na.rm = T) + N * mean(X_iso[, j], na.rm = T)^2 - 
                  sum(c(as.vector(n_array), N)) * mean.pool^2)/(sum(c(as.vector(n_array), 
                  N)) - 1))
              if (N == 1) 
                sd.pool <- sqrt((sum((as.vector(n_array) - 
                  1) * as.vector(SIG2_array[, j])) + as.vector(as.vector(n_array) %*% 
                  as.vector(MU_array[, j])^2) + N * mean(X_iso[, 
                  j], na.rm = T)^2 - sum(c(as.vector(n_array), 
                  N)) * mean.pool^2)/(sum(c(as.vector(n_array), 
                  N)) - 1))
              MU_array[, j] <- (MU_array[, j] - mean.pool)/sd.pool
              SIG2_array[, j] <- SIG2_array[, j]/sd.pool^2
          }
      }
      X_iso[, j] <- (X_iso[, j] - mean.pool)/sd.pool
      frac_mu[, j] <- frac_mu[, j]/sd.pool
      frac_sig2[, j] <- frac_sig2[, j]/sd.pool^2
  }
  R2jags::jags(jags.data, inits = jags.inits, parameters.to.save = jags.params, 
      model.file = model_filename, n.chains = mcmc$chains, 
      n.burnin = mcmc$burn, n.thin = mcmc$thin, n.iter = mcmc$chainLength, 
      DIC = mcmc$calcDIC)
}

make_jags_post_prop_long <- function(modfit, mu_tab, synth_df, n, ...) {
  ordered_sources <- sort(mu_tab$source)
  names(ordered_sources) <- as.character(seq_along(ordered_sources))
  dplyr::left_join(
    modfit$BUGSoutput$sims.matrix |>
      apply(2, ggdist::mean_hdci) |>
      dplyr::bind_rows(.id = "Variable") |>
      dplyr::filter(grepl("p\\.ind\\[", Variable)) |>
      tidyr::separate_wider_regex(
        Variable, c("p\\.ind\\[", N = "\\d+", "\\,", source = "\\d+", "\\]")
      ) |>
      dplyr::mutate(dplyr::across(c(N, source), as.numeric)) |>
      dplyr::mutate(source = .env$ordered_sources[source]) |>
      dplyr::select(source, N, Predicted = y, ymin, ymax),
    MixMustR::reshape_ref_data(synth_df, order_ref = mu_tab$source, ...) |>
      data.frame(check.names = FALSE) |>
      dplyr::select(-Unsampled) |>
      dplyr::mutate(N = seq_len(n())) |>
      tidyr::pivot_longer(
        cols = !N, names_to = "source", values_to = "Observed"
      ), by = join_by(N, source)
  ) |>
    dplyr::mutate(`Variant:` = as.character(.env$n))
}

make_jags_post_rand_long <- function(modfit, mu_tab, synth_df) {
  reff_df <- synth_df$sim_data |>
    dplyr::mutate(R = as.numeric(as.factor(paste0(Location, Habitat)))) |>
    dplyr::distinct(R, Location, Habitat) |>
    dplyr::arrange(R)
  dplyr::left_join(    
    modfit$BUGSoutput$sims.matrix |>
      colMeans() |>
      data.frame(mean = _) |>
      tibble::rownames_to_column("Variable") |>
      dplyr::filter(grepl("p\\.fac1\\[", Variable)) |>
      tidyr::separate_wider_regex(
        Variable, c("p\\.fac1\\[", R = "\\d+", "\\,", source = "\\d+", "\\]")
      ) |>
      dplyr::mutate(dplyr::across(c(R, source), as.numeric)) |>
      dplyr::select(mean, R, source) |>
      tidyr::pivot_wider(names_from = "source", values_from = "mean") |>
      MixMustR:::assign_new_names(c("R", mu_tab$Group)) |>
      dplyr::arrange(R) |>
      dplyr::left_join(reff_df, by = "R") |>
      tidyr::pivot_longer(
        cols = !c(R, Location, Habitat), names_to = "source",
        values_to = "Predicted"
      ),
    reshape_ref_data(synth_df, order_ref = mu_tab$Group) |>
      abs_log() |>
      data.frame(check.names = FALSE) |>
      dplyr::select(-Unsampled) |>
      (\(x, y)cbind(x, y))(
        x = dplyr::select(synth_df$props_df, Location, Habitat)
      ) |>
      dplyr::mutate(R = as.numeric(as.factor(paste0(Location, Habitat)))) |>
      tidyr::pivot_longer(
        cols = !c(R, Location, Habitat), names_to = "source",
        values_to = "Observed"
      ) |>
      dplyr::group_by(R, Location, Habitat, source) |>
      dplyr::summarise(Observed = exp(mean(Observed))) |>
      dplyr::ungroup(),
    by = join_by(R, Location, Habitat, source)
  )
}

reshape_jags_output <- function(model) {
  posts <- model$BUGSoutput$sims.matrix
  posts <- posts[, grepl("^p\\.ind\\[", colnames(posts))]
  matches <- regexec("p\\.ind\\[(\\d+),(\\d+)\\]", colnames(posts)) |>
    regmatches(colnames(posts), m = _)
  obs <- sapply(matches, "[[", 2) |> as.numeric()
  groups <- sapply(matches, "[[", 3) |> as.numeric()
  max_obs <- max(obs)
  max_groups <- max(groups)
  final_dims <- c(nrow(posts), max_obs, max_groups)
  final_array <- array(0, dim = final_dims)
  for (j in seq_along(unique(groups))) {
    final_array[, , j] <- posts[, groups == unique(groups)[j]]
  }
  final_array
}

mixsiar_bayes_R2 <- function(jagsfit, summary = TRUE, ...) {
  ref_data <- MixMustR::reshape_ref_data(...)
  ref_data <- ref_data[, colnames(ref_data) != "Unsampled"]
  all_posts <- reshape_jags_output(jagsfit)
  new_dims <- c(ncol(all_posts), dim(all_posts)[3], nrow(all_posts))
  pred_y <- res_y <- array(0, dim = new_dims)
  for (i in seq_len(dim(all_posts)[3])) {
    for (j in seq_len(nrow(all_posts))) {
      pred_y[, i, j] <- all_posts[j, , i]
    }
  }
  for (i in seq_len(dim(res_y)[3])) {
    res_y[, , i] <- ref_data - pred_y[, , i]
  }
  out_r2 <- matrix(0, dim(res_y)[3], ncol(res_y))
  colnames(out_r2) <- colnames(ref_data)
  for (i in seq_len(nrow(out_r2))) {
    pred_i <- pred_y[, , i]
    res_i <- res_y[, , i]
    var_pred <- apply(pred_i, 2, var)
    var_res <- apply(res_i, 2, var)
    out_r2[i, ] <- var_pred / (var_pred + var_res)
  }
  if (summary) {
    apply(out_r2, 2, ggdist::mean_hdci) |>
      dplyr::bind_rows(.id = "Source") |>
      dplyr::select(!(`.width`:`.interval`)) |>
      rename(mean = y, `2.5%HDI` = ymin, `97.5%HDI` = ymax) |>
      dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~round(.x, 2)))
  } else {
    out_r2
  }
}

extract_summary_r2 <- function(model_list, to_keep, ...) {
  purrr::map(model_list[to_keep], "model") |>
    purrr::map(function(x, ...) {
      wrangle_r2(x, ..., rm = FALSE)
    }, ...) |>
    dplyr::bind_rows() |>
    dplyr::mutate(Source = gsub("^Terrestrial", "T.", Source)) |>
    dplyr::rename(Mean = mean)
}

make_tbl_1 <- function(model_list_c, synth_list_c, model_list_d, synth_list_d,
                       jagsfits, ...) {
  to_filter <- which(
    MixMustR::mixmustr_models$sample_tracer &
      MixMustR::mixmustr_models$fix_unsampled &
      MixMustR::mixmustr_models$hierarchical
  )
  all_r2s_c <- extract_summary_r2(
    model_list = model_list_c, data_streams_list = synth_list_c,
    to_keep = to_filter, ...
  ) |>
    dplyr::mutate(Package = c(rep("", 3), "MixMustR", rep("", 3)),
                  Dataset = c(rep("", 3), "Convergent (1)", rep("", 3)))
  all_r2s_d <- extract_summary_r2(
    model_list = model_list_d, data_streams_list = synth_list_d,
    to_keep = to_filter, ...
  ) |>
    dplyr::mutate(Package = c(rep("", 3), "MixMustR", rep("", 3)),
                  Dataset = c(rep("", 3), "Divergent (2)", rep("", 3)))
  all_r2s_j <- mixsiar_bayes_R2(
    jagsfits$model, data_streams_list = synth_list_c, ...
  ) |>
    dplyr::mutate(
      Source = gsub("^Terrestrial", "T.", Source),
      Package = c(rep("", 2), "MixSIAR", rep("", 3)),
      Dataset = c(rep("", 2), "Convergent (1)", rep("", 3))
    ) |>
    dplyr::rename(Mean = mean)
  rbind(all_r2s_c, all_r2s_d, all_r2s_j) |>
    dplyr::relocate(Package, Dataset, .before = Source)
}

prior_uncertainty_plot <- function(x) {
  ln_sigma_x <- rep(0.1, length(x))
  plot_a <- MixMustR::evaluate_uncertainty(x, ln_sigma_x)$plot +
    labs(
      title = substitute(italic(pi) == "{" * a * "}", list(a = paste0(x, collapse = ", "))),
      subtitle = substitute(italic(sigma)["ln" * italic(rho)] == "{" * a * "}", list(a = paste0(ln_sigma_x, collapse = ", ")))
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
# SUPPLEMENTARY INFORMATION
# AND DIAGNOSTICS
###########################
###########################
###########################
elapsed_in_minutes <- function(model) {
  rstan::get_elapsed_time(model) |>
    rowSums() |>
    (\(x)round(x  / 60, 1))()
}

eval_sampler_diags <- function(model, ...) {
  treedpt <- function(x, target = 10) max(x[, "treedepth__"]) >= target
  divtrs <- function(x) mean(x[, "divergent__"]) * nrow(x)
  sampler_params <- rstan::get_sampler_params(model, inc_warmup = FALSE)
  n_mt <- sapply(sampler_params, treedpt, ...) |>
    sum()
  n_dv <- sapply(sampler_params, divtrs) |>
    sum() |>
    format(big.mark = ",")
  list(n_mt = n_mt, n_dv = n_dv)
}

plot_supp_checks <- function(model, n) {
  a <- plot_logit_probs_traces(model)
  b <- plot_logit_probs(model)
  c <- plot_rhat_and_eff(model)
  (a / b / c) +
    patchwork::plot_annotation(title = paste0("Model variant ", n))
}

plot_logit_probs_traces <- function(model) {
  to_string <- ggplot2::as_labeller(
    c(`logit(p[1,1])` = "Macroalgae", `logit(p[1,2])` = "Mangrove",
      `logit(p[1,3])` = "Saltmarsh", `logit(p[1,4])` = "Seagrass",
      `logit(p[1,5])` = "SPOM", `logit(p[1,6])` = "Terrestrial grass",
      `logit(p[1,7])` = "Unsampled")
  )
  bayesplot::mcmc_trace(
    model, pars = paste0("p[1,", 1:7, "]"), transformations = "logit",
    facet_args = list(nrow = 2, labeller = to_string)
  ) +
    theme_classic() +
    theme(legend.position = "bottom")
}

plot_logit_probs <- function(model) {
  to_string <- ggplot2::as_labeller(
    c(`logit(p[1,1])` = "Macroalgae", `logit(p[1,2])` = "Mangrove",
      `logit(p[1,3])` = "Saltmarsh", `logit(p[1,4])` = "Seagrass",
      `logit(p[1,5])` = "SPOM", `logit(p[1,6])` = "Terrestrial grass",
      `logit(p[1,7])` = "Unsampled")
  )
  bayesplot::mcmc_dens_overlay(
    model, pars = paste0("p[1,", 1:7, "]"), transformations = "logit",
    facet_args = list(nrow = 2, labeller = to_string)
  ) +
    theme_classic() +
    theme(legend.position = "bottom")
}

plot_rhat_and_eff <- function(model) {
  s_sum <- summary(model)$summary
  s_sum <- s_sum[grep("^p\\[", rownames(s_sum)), ]
  neffs <- data.frame(
    `\"# ESS\"` = s_sum[, "n_eff"],
    `widehat(italic(\"R\"))` = s_sum[, "Rhat"], check.names = FALSE
  ) |>
    tidyr::pivot_longer(`\"# ESS\"`:`widehat(italic(\"R\"))`,
                        names_to = "metric")
  vlines <- data.frame(x = 1, metric = "widehat(italic(\"R\"))")
  ggplot(data = neffs) +
    geom_histogram(mapping = aes(x = value), colour = "black",
                   position = "stack") +
    geom_vline(data = vlines, mapping = aes(xintercept = x),
               linetype = 2, colour = "grey30") +
    labs(x = "Metric") +
    facet_wrap(~metric, scales = "free", labeller = label_parsed) +
    theme_classic()
}

make_runtime_table <- function(model_list, tag) {
  cbind(
    MixMustR::mixmustr_models |>
      dplyr::select(!code_path) |>
      dplyr::mutate(across(everything(), ~if_else(.x, "Yes", "No"))) |>
      dplyr::rename(
        `Estimate tracer` = sample_tracer, `Fix unsampled` = fix_unsampled,
        `Hierarchical` = hierarchical
      ),
    purrr::map(model_list, "model") |>
      purrr::map(elapsed_in_minutes) |>
      dplyr::bind_rows() |>
      dplyr::rename_with(
        ~Hmisc::capitalize(gsub(":", " ", .x, fixed = TRUE))
      )
  ) |>
    dplyr::mutate(Dataset = c(tag, rep("", n() - 1))) |>
    dplyr::relocate(Dataset, .before = `Estimate tracer`)
}

make_table_s2 <- function(model_list_c, model_list_d) {
  rbind(make_runtime_table(model_list_c, "Convergent"),
        make_runtime_table(model_list_d, "Divergent"))
}

make_fig_s1_8 <- function(to_filter, model_list, ylab, subt, ...) {
  rbind(
    lapply(to_filter, function(x, mod_list, ...) {
      MixMustR::make_post_prop_long(
        mod_list[[x]]$model, ..., n = x
      )
    }, mod_list = model_list, ...) |>
      dplyr::bind_rows()
  ) |>
    MixMustR::plot_multiple_faceted_scatter_avg() +
    labs(y = ylab,
         title = "Model fit performance, mixing proportions",
         subtitle = subt) +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 8), 
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12)
          ) +
    guides(shape = guide_legend(override.aes = list(size = 3.5)),
           fill = guide_legend(override.aes = list(alpha = 1)))
}

run_prior_model_si <- function(data_streams_list, tracer_pars, model_path) {
  sdata <- MixMustR::mixmustr_wrangle_input(
    data_streams_list, tracer_pars, model_path, sigma_ln_rho = 1,
    sample_tracer = FALSE, fix_unsampled = TRUE, hierarchical = FALSE
  )
  mod_name_suffix_i <- tools::file_path_sans_ext(basename(model_path))
  timing <- system.time({
    model <- MixMustR::run_mixmod(
      model_path, mod_name_suffix_i, data = sdata, iter = 1e4, warmup = 5e3,
      chains = 4, cores = 4
    )
  })
  models <- vector(mode = "list", length = 6)
  models[[6]] <- list(timing = timing, model = model)
  names(models)[6] <- mod_name_suffix_i
  models
}
