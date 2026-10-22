MixMustR description <img src="docs/logo.png" width = 180 alt="MixMustR logo" align="right" />
=========================================================================================

[![Ask Us Anything][0a]][0b]
![Open Source Love][0c]
[![license](https://img.shields.io/badge/license-MIT%20+%20file%20LICENSE-lightgrey.svg)](https://choosealicense.com/)

[0a]: https://img.shields.io/badge/Ask%20us-anything-1abc9c.svg
[0b]: https://github.com/open-aims/MixMustR_paper/issues/new
[0c]: https://badges.frapsoft.com/os/v2/open-source.svg?v=103

# MixMustR paper

# Package

https://github.com/open-aims/MixMustR

## Bayesian mixture model

## Instructions

All analyses were done in `R`. This routine uses the [targets](https://github.com/ropensci/targets) R package to produce all outputs. First install `targets`:

```r
install.packages("targets")
```

Next you need to open an R session with working directory set to the root of the project.

This routine uses some packages detailed under `_targets.R` (plus `gt`, used in `docs/Supps.qmd`), so make sure to install them before running anything. Rendering the manuscript (`docs/MS.qmd`) and supplementary information (`docs/Supps.qmd`) also requires [Quarto](https://quarto.org/) and a LaTeX distribution with `xelatex`.

Then, to generate the analyses and output visualisations simply run:

```r
source("run.R")
```

Then, to load in each of the targets from `_targets.R` simply run in an R session:

```
targets::tar_load(target_name)
```

## HPC simulation study (read-only)

The simulation study summarised in Figs 4–5 and Figs S1–S7 was run on a SLURM HPC cluster. The code is provided under `hpc/` (plus `R/functions_hpc.R`) for transparency only; it is **not** part of the `targets` pipeline and is not meant to be executed from this repository:

- `hpc/Dockerfile`: container image used on the cluster.
- `hpc/make_hpc_pack_and_jobs.R`: builds the simulation inputs and the job table `hpc/jobs_dev.csv`.
- `hpc/run_csv_array_dev.sbatch`: SLURM array script (one task per row of `hpc/jobs_dev.csv`).
- `hpc/worker_dev.R`: per-task worker that simulates data, fits a MixMustR model and saves diagnostics.

The HPC outputs required by the `targets` pipeline are included so that the manuscript and supplementary information can be reproduced without access to an HPC:

- `hpc/jobs_dev.csv`: simulation design.
- `hpc/diagnostics_rds/diag_model_*.rds`: per-fit diagnostics (6,400 files), read by the `hpc_results` target.

## License

This repository is provided by the authors under the MIT License ([MIT](http://opensource.org/licenses/MIT)).

## Bug reporting
* Please [report any issues or bugs](https://github.com/open-aims/MixMustR_paper/issues).
