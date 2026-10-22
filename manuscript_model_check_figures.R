# Minimal script for the three model-check figures in the manuscript:
#
#   dispersion_by_gene.pdf                   Pearson dispersion ratio per gene
#   interval_check_by_gene_concentration.pdf observed mean and residual SD
#                                            against the 95% simulated interval,
#                                            per gene and concentration
#   check_fit_by_gene.pdf                    the same intervals, one page per gene
#
# This is the subset of model_checks_brms_multivariate.R that these figures
# depend on, with the same seed and the same calculations, so it reproduces
# the figures in model_checks/. The rationale for each calculation is given in
# full in that script.
#
# Inputs are the two multivariate brms fits (78 genes, cubic polynomial in
# log concentration, 4 chains of 1000 post-warmup draws):
#   fit3k_bb_DRomicsGenes_poly_iter1000.RData     beta-binomial
#   fit3k_binom_DRomicsGenes_poly_iter1000.RData  binomial
# The beta-binomial fit is produced in fit_brms_betaBinomial_poly.R.
#
# Runs in a few minutes on a desktop; posterior_predict() over all 78 genes is
# the slowest step.

library(brms)
library(tidyverse)

# The seed is set again immediately before each fit's posterior_predict() call
# in fit_arrays(), so a fit's simulated data do not depend on which fits were
# processed before it.
set.seed(10)
check_seed <- 10

out_dir <- "manuscript_figures"
dir.create(out_dir, showWarnings = FALSE)

# The beta-binomial is first: genes are ordered by its dispersion ratio, and
# it takes the first colour.
fit_files <- c(
  beta_binomial = "fit3k_bb_DRomicsGenes_poly_iter1000.RData",
  binomial = "fit3k_binom_DRomicsGenes_poly_iter1000.RData"
)

####---- extracting what the checks need from a fit ----####

load_fit <- function(file) {
  e <- new.env()
  nm <- load(file, envir = e)
  if (length(nm) != 1 || !inherits(e[[nm]], "brmsfit")) {
    stop(file, " should hold exactly one brmsfit.", call. = FALSE)
  }
  e[[nm]]
}

# brms strips the underscore from a response variable when naming it, so
# "gene_Amillepora00288" is the response "geneAmillepora00288".
gene_label <- function(resp) sub("^gene_?", "", resp)

# Variance of the observed count given the fitted mean count m. brms returns
# the expected count from posterior_epred() for the binomial families, so the
# probability is m / trials. The beta-binomial variance is
# n p (1 - p) (phi + n) / (phi + 1) in the mean-precision form brms uses.
obs_variance <- function(family, m, trials, phi) {
  p <- m / trials
  switch(
    family,
    binomial = trials * p * (1 - p),
    beta_binomial = trials * p * (1 - p) * (phi + trials) / (phi + 1),
    stop("No variance is defined for family '", family, "'.", call. = FALSE)
  )
}

# Everything the checks need, one element per response. Predictions are drawn
# for the whole multivariate fit in one call rather than once per response,
# which would rebuild the predictions for all 78 genes each time.
fit_arrays <- function(fit, seed) {
  resps <- fit$formula$responses
  sdat <- brms::standata(fit)
  draws <- posterior::as_draws_matrix(fit)
  set.seed(seed)
  yrep <- brms::posterior_predict(fit)
  epred <- brms::posterior_epred(fit)
  out <- lapply(seq_along(resps), function(i) {
    r <- resps[i]
    fam <- fit$family[[r]]$family
    y <- sdat[[paste0("Y_", r)]]
    # phi is a single parameter per gene in these fits, so its draws are one
    # column, recycled across observations.
    phi <- if (fam == "beta_binomial") {
      matrix(draws[, paste0("phi_", r)], nrow(draws), length(y))
    }
    list(
      resp = r,
      family = fam,
      y = y,
      trials = sdat[[paste0("trials_", r)]],
      phi = phi,
      yrep = yrep[,, i],
      epred = epred[,, i]
    )
  })
  names(out) <- resps
  out
}

####---- the checks, per response ----####

# As bayesnec::dispersion(): per posterior draw, the Pearson residual sum of
# squares of the observed data divided by that of data simulated from the same
# draw. A ratio near 1 is what the fitted model expects. Terms with a
# non-positive fitted variance are excluded rather than allowed to give NaN;
# none arise in these fits.
dispersion_check <- function(a) {
  n_draw <- nrow(a$epred)
  trials <- matrix(a$trials, n_draw, ncol(a$epred), byrow = TRUE)
  v <- obs_variance(a$family, a$epred, trials, a$phi)
  y_mat <- matrix(a$y, n_draw, ncol(a$epred), byrow = TRUE)
  v[!is.finite(v) | v <= 0] <- NA_real_
  obs_ss <- rowSums((y_mat - a$epred)^2 / v, na.rm = TRUE)
  sim_ss <- rowSums((a$yrep - a$epred)^2 / v, na.rm = TRUE)
  ratio <- obs_ss / sim_ss
  tibble::tibble(
    disp_median = stats::median(ratio),
    disp_lower_95 = unname(stats::quantile(ratio, 0.025)),
    disp_upper_95 = unname(stats::quantile(ratio, 0.975))
  )
}

# As bayesnec::check_fit_stats(): per concentration, the observed mean and
# residual SD, and the 95% interval of each over the simulated data sets. The
# residual SD is taken about the posterior median fitted mean, and the same
# fitted mean is used for the observed and the simulated residuals; using each
# draw's own mean for the simulated ones would remove the discrepancy being
# measured. In these fits every observation in a group shares one
# concentration and one number of trials, so the residual SD of a group equals
# the SD of its raw counts.
check_fit_groups <- function(a, log_x) {
  mu <- apply(a$epred, 2, stats::median)
  conc <- sort(unique(log_x))
  purrr::map_dfr(conc, function(g) {
    idx <- which(log_x == g)
    obs_y <- a$y[idx]
    sim_y <- a$yrep[, idx, drop = FALSE]
    mean_draws <- rowMeans(sim_y)
    sd_draws <- apply(sweep(sim_y, 2, mu[idx], "-"), 1, stats::sd)
    tibble::tibble(
      log_x = g,
      obs_mean = mean(obs_y),
      sim_mean_lo = unname(stats::quantile(mean_draws, 0.025)),
      sim_mean_hi = unname(stats::quantile(mean_draws, 0.975)),
      obs_sd = stats::sd(obs_y - mu[idx]),
      sim_sd_lo = unname(stats::quantile(sd_draws, 0.025)),
      sim_sd_hi = unname(stats::quantile(sd_draws, 0.975))
    )
  })
}

run_checks <- function(fit_name, file) {
  message("Checking ", fit_name, " (", file, ")")
  fit <- load_fit(file)
  log_x <- fit$data$log_x
  arrays <- fit_arrays(fit, check_seed)
  stack <- function(f) {
    purrr::imap_dfr(arrays, function(a, r) {
      dplyr::mutate(f(a), resp = r, .before = 1)
    })
  }
  list(
    dispersion = stack(dispersion_check),
    check_fit = stack(\(a) check_fit_groups(a, log_x))
  )
}

checks <- purrr::imap(fit_files, function(file, nm) {
  out <- run_checks(nm, file)
  gc()
  out
})

bind_part <- function(part) {
  purrr::imap_dfr(checks, function(x, nm) {
    dplyr::mutate(x[[part]], fit = nm, .before = 1)
  }) |>
    dplyr::mutate(gene = gene_label(resp), .after = fit)
}

dispersion_tab <- bind_part("dispersion")
check_fit_tab <- bind_part("check_fit") |>
  dplyr::mutate(conc = round(exp(log_x), 1))

####---- figures ----####

# Fixed colour per fit, with shape as a second encoding so the fits are not
# distinguished by colour alone.
fit_cols <- c(beta_binomial = "#2a78d6", binomial = "#eb6834")
fit_shapes <- c(beta_binomial = 16, binomial = 17)

# Genes ordered by the beta-binomial dispersion ratio, in all three figures.
gene_order <- dispersion_tab |>
  dplyr::filter(fit == names(fit_files)[1]) |>
  dplyr::arrange(disp_median) |>
  dplyr::pull(gene)

p_disp <- dispersion_tab |>
  dplyr::mutate(gene = factor(gene, levels = gene_order)) |>
  ggplot2::ggplot(ggplot2::aes(
    x = gene,
    y = disp_median,
    colour = fit,
    shape = fit
  )) +
  ggplot2::geom_hline(yintercept = 1, colour = "grey50", linewidth = 0.4) +
  ggplot2::geom_linerange(
    ggplot2::aes(ymin = disp_lower_95, ymax = disp_upper_95),
    linewidth = 0.5
  ) +
  ggplot2::geom_point(size = 1.6) +
  ggplot2::scale_y_log10() +
  ggplot2::scale_colour_manual(values = fit_cols) +
  ggplot2::scale_shape_manual(values = fit_shapes) +
  ggplot2::labs(
    x = NULL,
    y = "Pearson dispersion ratio (observed / simulated)",
    colour = NULL,
    shape = NULL,
    caption = "Median and 95% interval over posterior draws; 1 is the value expected under the fitted model."
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 5
    ),
    legend.position = "top"
  )
ggplot2::ggsave(
  file.path(out_dir, "dispersion_by_gene.pdf"),
  p_disp,
  width = 25,
  height = 12,
  units = "cm"
)

# Each observed statistic classed as above, inside or below the 95% interval
# of the same statistic simulated from the fit. Outside the 95% interval
# corresponds to a ppp outside [0.025, 0.975], so about 5% of cells are
# expected outside under a correct model. Blue where the observed value is
# above the interval, red where it is below.
interval_levels <- c("above interval", "inside", "below interval")
interval_cols <- c(
  "above interval" = "#1c5cab",
  "inside" = "#e3e2de",
  "below interval" = "#c0392b"
)

# Long form, one row per fit, gene, concentration and statistic. Used by both
# the heatmap and the per-gene pages.
interval_tab <- dplyr::bind_rows(
  check_fit_tab |>
    dplyr::transmute(
      fit,
      gene,
      conc,
      statistic = "mean",
      observed = obs_mean,
      lo = sim_mean_lo,
      hi = sim_mean_hi
    ),
  check_fit_tab |>
    dplyr::transmute(
      fit,
      gene,
      conc,
      statistic = "residual SD",
      observed = obs_sd,
      lo = sim_sd_lo,
      hi = sim_sd_hi
    )
) |>
  dplyr::mutate(
    position = dplyr::case_when(
      observed > hi ~ "above interval",
      observed < lo ~ "below interval",
      .default = "inside"
    ),
    position = factor(position, levels = interval_levels)
  )

p_interval <- interval_tab |>
  dplyr::mutate(
    gene = factor(gene, levels = gene_order),
    conc = factor(conc)
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = conc, y = gene, fill = position)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
  ggplot2::scale_fill_manual(values = interval_cols, drop = FALSE) +
  ggplot2::facet_grid(cols = ggplot2::vars(statistic, fit)) +
  ggplot2::labs(
    x = "1-MN (ug/L)",
    y = NULL,
    fill = "observed value",
    caption = paste(
      "Observed statistic against the 95% interval of the same statistic",
      "simulated from the fit. About 5% of cells are expected outside under a correct model."
    )
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 5),
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "top"
  )
ggplot2::ggsave(
  file.path(out_dir, "interval_check_by_gene_concentration.pdf"),
  p_interval,
  width = 30,
  height = 26,
  units = "cm"
)

# One page per gene, the fits side by side, in the layout of
# bayesnec's plot.checkfit(). A point outside its bar is a concentration the
# fit does not reproduce.
plot_gene_check <- function(g) {
  dodge <- ggplot2::position_dodge(width = 0.5)
  interval_tab |>
    dplyr::filter(gene == g) |>
    dplyr::mutate(conc = factor(conc)) |>
    ggplot2::ggplot(ggplot2::aes(x = conc, colour = fit, shape = fit)) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = lo, ymax = hi),
      position = dodge,
      linewidth = 0.6
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = observed),
      position = dodge,
      size = 2
    ) +
    ggplot2::facet_wrap(ggplot2::vars(statistic), scales = "free_y") +
    ggplot2::scale_colour_manual(values = fit_cols) +
    ggplot2::scale_shape_manual(values = fit_shapes) +
    ggplot2::labs(
      title = g,
      x = "1-MN (ug/L)",
      y = "observed (point) and simulated 95% (bar)",
      colour = NULL,
      shape = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top")
}

grDevices::pdf(
  file.path(out_dir, "check_fit_by_gene.pdf"),
  width = 9,
  height = 4.5
)
for (g in gene_order) {
  print(plot_gene_check(g))
}
grDevices::dev.off()
