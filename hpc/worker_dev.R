library(MASS)
library(MixMustR)
library(tidyverse)
library(gtools)
library(rlang)
library(patchwork)
library(tools)
library(posterior)
library(future)
library(furrr)
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, args) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("Missing value for ", flag)
  args[[i + 1]]
}

sim_id <- as.integer(get_arg("--sim_id", args))
model_id <- as.integer(get_arg("--model_id", args))
truth_stream <- as.integer(get_arg("--truth_stream", args))
prop_seed <- as.integer(get_arg("--prop_seed", args))
tracer_seed <- as.integer(get_arg("--tracer_seed", args))
stan_seed <- as.integer(get_arg("--stan_seed", args))
input_file <- get_arg("--input_file", args)
M <- as.integer(get_arg("--M", args))
G <- as.integer(get_arg("--G", args))
degree <- get_arg("--degree", args)
noises <- get_arg("--noises", args) |>
  (\(x)strsplit(x, ";")[[1]])() |>
  as.numeric()
tracers <- get_arg("--tracers", args) |>
  (\(x)strsplit(x, ";")[[1]])()
source_names <- get_arg("--source_names", args) |>
  (\(x)strsplit(x, ";")[[1]])()
J <- as.integer(length(source_names))
model_chars <- MixMustR::mixmustr_models[model_id, ]

cat("Simulation ID :" , sim_id, "\n")
cat("Model ID :" , model_id, "\n")
cat("  Sample tracer :" , model_chars$sample_tracer, "\n")
cat("  Fix unsampled :" , model_chars$fix_unsampled, "\n")
cat("  Hierarchical :" , model_chars$hierarchical, "\n")
cat("Truth stream :" , truth_stream, "\n")
cat("Proportions seed :" , prop_seed, "\n")
cat("Tracer seed :" , tracer_seed, "\n")
cat("Stan seed :" , stan_seed, "\n")
cat("Input file :" , input_file, "\n")
cat("M :" , M, "\n")
cat("G :" , G, "\n")
cat("J :" , J, " (", paste0(source_names, collapse = ", "), ")", "\n")
cat("Tracers :" , paste0(tracers, collapse = ", "), "\n")
cat("Degree :" , degree, " (", paste0(noises, collapse = ", "), ")", "\n")

load(input_file)

for (i in seq_along(tracer_parameters)) {
  tracer_parameters[[i]] <- (tracer_parameters[[i]]) |>
    dplyr::filter(source %in% .env$source_names) |>
    dplyr::select(dplyr::all_of(c("source", tracers)))
}
bcs_cols <- c(
  "source", "Study",
  sapply(tracers, function(x) c(x, paste(x, "(SD)"), paste(x, "(n)")))
)
bcs <- bcs |>
  dplyr::filter(source %in% .env$source_names) |>
  dplyr::select(dplyr::all_of(intersect(names(bcs), bcs_cols)))
safe_sim_fct <- purrr::safely(sim_fct, otherwise = NULL)
res <- safe_sim_fct(
  prop_seed, tracer_seed, stan_seed, tracer_parameters, bcs, sim_id,
  model_id, truth_stream, to_write = TRUE, source_names = source_names,
  noises = noises, M = M, G = G, J = J
)
if (!is.null(res$error)) {
  err_tbl <- tibble::tibble(
    sim_id = sim_id, model_id = model_id,
    sample_tracer = model_chars$sample_tracer,
    fix_unsampled = model_chars$fix_unsampled,
    hierarchical = model_chars$hierarchical, truth_stream = truth_stream,
    prop_seed = prop_seed, tracer_seed = tracer_seed,
    stan_seed = stan_seed, M = M, G = G, J = J,
    source_names = paste0(source_names, collapse = ";"),
    noises = paste0(noises, collapse = ";"),
    error_message = conditionMessage(res$error)
  )
  out_name <- paste0("revision_sim_", sim_id, "_ERROR.rds")
  saveRDS(err_tbl, file = file.path("copilot", "sim_models", out_name))
  cat(sim_id, " has failed.")
} else {
  out <- purrr::map(
    res$result,
    ~ dplyr::mutate(.x, sim_id = sim_id, .before = 1)
  )
  out_name <- paste0("revision_sim_", sim_id, ".rds")
  saveRDS(out, file = file.path("copilot", "sim_models", out_name))
  message("Successfully finished simulation ", sim_id, ".")
}
