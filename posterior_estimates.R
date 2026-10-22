# Posterior NSEC estimates for the 78 DRomics-selected genes
#
# Fits a multivariate beta-binomial cubic polynomial (brms) to the
# batch-corrected counts of the 78 genes retained by DRomics (FDR 0.05),
# estimates the no-significant-effect concentration (NSEC) of each gene in
# both directions with toxval::nsec_multi(), and adds the posterior median of
# the expected count in the control. Writes posterior_estimates.csv.
#
# Inputs (in the working directory):
#   BECdata_wSampleID.csv          batch-corrected counts, one row per gene
#   DromicsFittedGenesFDR0.05.csv  GeneID of the 78 genes
#   Metadata_1-MN.csv              sample ID and measured 1-MN (ug/L)
#
# Software. The published fit was run with brms 2.20.4 and rstan 2.32.3.
# toxval is installed from GitHub at the commit used for the manuscript:
#   remotes::install_github(
#     "open-AIMS/toxval@41b7f789fedda1982743e17f6701fd9aa926633f"
#   )
# This is the `pre-refactor` tag; later versions of toxval change the
# nsec_multi() interface.
#
# Reproducibility. Stan is given seed = 100, but identical draws are only
# obtained with the same Stan, rstan and brms versions, compiler and platform.
# brms 2.23 writes the same model with a vectorised likelihood, which changes
# the Stan program and therefore the draws. On another set-up the estimates
# differ by Monte Carlo error. To reproduce the
# published values exactly, place the saved fit (fit_file below) in the working
# directory and the script loads it instead of refitting.

library(brms)
library(toxval)
library(tidyverse)

fit_file <- "fit3k_bb_DRomicsGenes_poly_iter1000.RData"

####---- data ----####

bec_dat <- read.csv("BECdata_wSampleID.csv")

# Total count of each gene over the 48 samples, used as the number of trials,
# so the beta-binomial models each sample's share of that gene's total.
gene_counts <- read.csv("DromicsFittedGenesFDR0.05.csv") |>
  dplyr::inner_join(bec_dat, by = "GeneID") |>
  dplyr::mutate(max = rowSums(dplyr::across(dplyr::where(is.numeric)))) |>
  tidyr::pivot_longer(
    cols = dplyr::starts_with("X"),
    names_to = "SampleID",
    values_to = "y"
  )

# The controls have a measured concentration of 0, which is set to 0.3 ug/L
# so that it can be log-transformed.
dat <- read.csv("Metadata_1-MN.csv") |>
  dplyr::left_join(gene_counts, by = "SampleID") |>
  dplyr::mutate(
    x = as.numeric(as.character(Measured)),
    x = ifelse(x == 0, 0.3, x),
    log_x = log(x)
  )

dat_trials <- dat |>
  dplyr::select(-y) |>
  tidyr::pivot_wider(
    names_from = GeneID,
    values_from = max,
    names_prefix = "n_"
  )

dat_use <- dat |>
  dplyr::select(-max) |>
  tidyr::pivot_wider(
    names_from = GeneID,
    values_from = y,
    names_prefix = "gene_"
  ) |>
  dplyr::left_join(
    dat_trials,
    by = c("SampleID", "Tr", "Settlement", "Nominal", "Measured", "x", "log_x")
  )

gene_ids <- unique(dat$GeneID)

####---- model ----####

# One response per gene, each with its own cubic in log concentration and its
# own precision; residual correlations are not modelled.
bf_3k <- purrr::map(gene_ids, \(g) {
  brms::bf(stats::as.formula(paste0(
    "gene_", g, " | trials(n_", g, ") ~ poly(log_x, 3, raw = TRUE)"
  )))
}) |>
  purrr::reduce(`+`)

if (file.exists(fit_file)) {
  load(fit_file)
  # Guards against pairing a saved fit with different input data.
  fit_cols <- setdiff(colnames(fit3k$data), "poly(log_x, 3, raw = TRUE)")
  stopifnot(isTRUE(all.equal(
    as.data.frame(fit3k$data[, fit_cols]),
    as.data.frame(dat_use[, fit_cols]),
    check.attributes = FALSE
  )))
} else {
  iter <- 2000
  fit3k <- brms::brm(
    bf_3k + brms::set_rescor(FALSE),
    data = dat_use,
    family = "beta_binomial",
    chains = 4,
    cores = 4,
    init = 0,
    seed = 100,
    iter = iter,
    warmup = iter - 1000
  )
  save(fit3k, file = fit_file)
}

####---- NSEC ----####

# NSEC references are the 1st (decreasing) and 99th (increasing) percentiles
# of the control posterior, evaluated on a 20-point grid of log concentration.
# A draw that never crosses its reference is assigned the highest
# concentration, 882 ug/L, so an estimate of 882 means no effect was detected
# within the tested range.
nsec_both <- toxval::nsec_multi(
  object = fit3k,
  resolution = 20,
  type = "both",
  posterior = FALSE,
  xform = exp,
  x_var = "log_x",
  trials_var = "n_",
  multi_var = "gene_"
)

####---- control ----####

# The control is predicted at the lowest x with each gene's trials, as
# nsec_multi() does for the first point of its grid, so ctrl_val is on the
# same scale as the references dec.1% and inc.99%.
trials_data <- fit3k$data |>
  dplyr::select(dplyr::starts_with("n_")) |>
  unique()
pred_dat <- data.frame(log_x = min(fit3k$data$log_x)) |>
  cbind(trials_data, row.names = NULL)

control_draws <- brms::posterior_epred(
  fit3k,
  newdata = pred_dat,
  re_formula = NA
)[, 1, ]

# brms drops underscores from response names, so genes are matched on the
# stripped name rather than by position.
gene_vars <- paste0("gene_", gene_ids)
colnames(control_draws) <- gene_vars[
  match(colnames(control_draws), gsub("_", "", gene_vars))
]

control <- tibble::tibble(
  vars = colnames(control_draws),
  ctrl_val = apply(control_draws, 2, stats::median)
)

posterior_estimates <- nsec_both |>
  dplyr::left_join(control, by = "vars")

stopifnot(
  nrow(posterior_estimates) == 78,
  !anyNA(posterior_estimates$ctrl_val)
)

write.csv(
  posterior_estimates,
  file = "posterior_estimates.csv",
  row.names = FALSE
)
