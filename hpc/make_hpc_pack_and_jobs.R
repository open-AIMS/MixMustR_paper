library(MASS)
library(MixMustR)
library(tidyverse)
library(gtools)
library(rlang)
library(patchwork)
library(tools)
library(MixSIAR)
library(posterior)
library(future)
library(furrr)
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("R/functions_hpc.R")
#############################
#############################
#############################
# READ IN GLOBAL PARAMETERS #
#############################
#############################
#############################
synthetic_data <- data.frame(group = paste0("G", rep(1:10, each = 10)))
tracer_parameters <- MixMustR::tracer_parameters
# Simulate a outside-of-convex hill tracer value for the unsampled
S <- tracer_parameters$mus |>
  tibble::column_to_rownames("source") |>
  as.matrix()
S_scaled <- scale(S)
center_tracer <- attr(S_scaled, "scaled:center")
scale_tracer  <- attr(S_scaled, "scaled:scale")
# Direction of strongest source separation in scaled tracer space
pc <- prcomp(S_scaled, center = FALSE, scale. = FALSE)
# Function to create an off-centroid unsampled signature
# unsampled_centroid <- make_unsampled_signature(delta = 0)
# unsampled_near     <- make_unsampled_signature(delta = 0.75)
# unsampled_mid      <- make_unsampled_signature(delta = 1.5)
unsampled_far      <- make_unsampled_signature(delta = 3)
# `Study` here is meaningless, just kept so `MixMustR:::reshape_isotope_df`
# doesn't fail
bcs <- dplyr::left_join(
  MixMustR::bcs_si, dplyr::select(MixMustR::bcs_fa, !c(Study, Taxa)),
  by = "source"
) |>
  dplyr::mutate(
    dplyr::across(contains("SD"), ~1), dplyr::across(contains("(n)"), ~10L)
  )
# Simulate a median SD value for the unsampled
# normalise using log so as not to over-inflate the simulation
tracers <- names(dplyr::select(tracer_parameters$mus, !source))
unsmpd <- tibble::tibble(source = "Unsampled", Study = "Simulated")
for (i in seq_along(tracers)) {
  trc <- tracers[i]
  tr_sd <- sqrt(exp(mean(log(bcs[[paste0(trc, " (SD)")]]^2), na.rm = TRUE)))
  unsmpd[[paste0(trc, " (SD)")]] <- tr_sd
  unsmpd[[paste0(trc, " (n)")]] <- as.integer(round(mean(bcs[[paste0(trc, " (n)")]], na.rm = TRUE)))
}
bcs <- dplyr::left_join(unsampled_far, unsmpd, by = "source") |>
  dplyr::bind_rows(bcs, ... = _)
tracer_parameters$sigmas <- tracer_parameters$sigmas |>
  dplyr::mutate(dplyr::across(!matches("source"), ~1))
tracer_parameters$ns <- tracer_parameters$ns |>
  dplyr::mutate(dplyr::across(!matches("source"), ~10L))
rm(tracers, trc, tr_sd, unsmpd)
save.image("hpc/hpc_pack.RData")

# simulate many proportion datasets and generation processes
pnc <- function(...) paste0(..., collapse = ";")
three_sources <- c("Mangrove", "SPOM", "Terrestrial grass", "Unsampled") |>
  pnc()
seven_sources <- c("Macroalgae", "Mangrove", "Saltmarsh", "Seagrass",
                   "SPOM", "Terrestrial grass", "Unsampled") |>
  pnc()
two_tracers <- c("d(13C/12C)", "d(15N/14N)") |>
  pnc()
all_tracers <- (names(tracer_parameters$mus)[-1]) |>
  pnc()
n_sims <- 50
expand.grid(
  tracers = c(two_tracers, all_tracers),
  source_names = c(three_sources, seven_sources), degree = c("con", "div"),
  sim_id = seq_len(n_sims), model_id = 1:8, truth_stream = 1:2, G = 10, M = 10,
  input_file = "hpc_pack.RData"
) |>
  dplyr::mutate(
    sim_id = seq_len(dplyr::n()), prop_seed = 10000L + sim_id,
    tracer_seed = 20000L + sim_id, stan_seed = 30000L + sim_id,
    noises = dplyr::case_when(
      degree == "con" ~ "0.5;1", degree == "div" ~ "0.5;0.001",
      .default = NA
    )
  ) |>
  write.csv("hpc/jobs_dev.csv", row.names = FALSE, quote = FALSE)
