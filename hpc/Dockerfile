FROM rocker/r-ver:4.5.3

## Install packages
RUN apt-get update \ 
  && apt-get install -y --no-install-recommends \ 
    libgdal-dev \ 
    libproj-dev \ 
    libgeos-dev \ 
    lbzip2 \ 
    libnetcdf-dev \ 
    libsqlite3-dev \ 
    libssl-dev \ 
    libudunits2-dev \ 
    netcdf-bin \ 
    sqlite3 \ 
    imagemagick \ 
    make \ 
    ghostscript \ 
    poppler-utils \ 
    zip \ 
    wget \ 
    fonts-dejavu-extra \ 
    curl \ 
    tk \ 
    openjdk-11-jre \ 
    gdebi-core \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    cmake \
  && rm -rf /var/lib/apt/lists/* 

## Install R packages from PPPM (based on a date - YYYY-MM-DD)
RUN R -e "options(repos = \
    list(CRAN = 'https://packagemanager.posit.co/cran/2026-03-01/')); \
  all_pkgs <- installed.packages()[, 'Package', drop = TRUE]; \
  to_install <- c('remotes', 'tidyverse', 'rlang', 'rstan', 'abind', 'beepr', 'brms', 'data.table', 'digest', 'JuliaCall', 'magrittr', 'matrixStats', 'plotly', 'popbio', 'R.matlab', 'readxl', 'tictoc', 'dplyr', 'usethis', 'knitr', 'rmarkdown', 'testthat', 'sf', 'ggdist', 'patchwork', 'cowplot', 'gtools', 'rjags', 'R2jags', 'MixSIAR', 'furrr'); \
  for (i in seq_along(to_install)) {; \
    if (!to_install[i] %in% all_pkgs) {; \
      install.packages(to_install[i]); \
      all_pkgs <- installed.packages()[, 'Package', drop = TRUE]; \
    }; \
  }; \
  install.packages('cmdstanr', repos = c('https://mc-stan.org/r-packages/', getOption('repos'))); \
  library(cmdstanr); \
  check_cmdstan_toolchain(); \
  install_cmdstan(cores = 2); \
  remotes::install_github('open-aims/MixMustR', ref = 'dev'); \
"

WORKDIR /home/Project
