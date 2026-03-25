rm(list = ls(all = TRUE))  ## efface les données
source('~/thib/projects/tools/R_lib.r')
source('utils.r')
setwd('/Users/thibault/thib/projects/reliable_info/reliable_info_identif/')

library(tidyverse)
library(posterior)
library(bayesplot)
library(loo)
library(knitr)

exp <- 'exp15'
models <- c('linear_report', 'BLO_report', 'BLO_noise_report')

# -------------------------
# LOAD DATA
# -------------------------
load(paste0("./data/data_list_proba_", exp, ".rdata"))
if (exists("data_list_report")) {
  data_list <- data_list_report
} else if (exists("data_list")) {
  data_list <- data_list
} else {
  stop("No `data_list_report` or `data_list` found in the loaded .rdata.")
}
N <- data_list$N

# ---------------------------------------------------
## LOAD MODELS FITS  (needed BEFORE diagnosis)
# ---------------------------------------------------
fits <- vector("list", length(models))
for (i in seq_along(models)) {
  load(paste0('./results/fits/', exp, '/fit_', models[i], '_', exp, '.rdata'))
  fits[[i]] <- fit
}

dir.create(paste0("./results/plots/", exp), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0("./results/summary/", exp), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------
## MODEL DIAGNOSIS
# ---------------------------------------------------

## Pairs plots (mu_* + mu_pr[*])
for (i in seq_along(models)) {
  fit <- fits[[i]]
  posterior_df <- as_draws_df(fit$draws())
  
   selected_params <- posterior_df[, grepl("^mu", colnames(posterior_df)) & !grepl("^mu_pr", colnames(posterior_df))]
  ##keep <- grepl('^mu_', colnames(posterior_df)) | grepl('^mu_pr\\[', colnames(posterior_df))
  ##selected_params <- posterior_df[, keep, drop = FALSE]

  plot_pairs <- mcmc_pairs(selected_params)
  plot_file <- paste0("./results/plots/", exp, "/pairs_plot_", models[i], "_", exp, ".pdf")
  ggsave(filename = plot_file, plot = plot_pairs, width = 7, height = 7)
}

## Traces (mu_* + mu_pr[*])
for (i in seq_along(models)) {
  fit <- fits[[i]]
  posterior_df <- as_draws_df(fit$draws())

  keep <- grepl('^mu_', colnames(posterior_df)) | grepl('^mu_pr\\[', colnames(posterior_df))
  selected_params <- posterior_df[, keep, drop = FALSE]

  plot_trace <- mcmc_trace(selected_params)
  plot_file <- paste0('./results/plots/', exp, '/trace_plot_', models[i], "_", exp, ".pdf")
  ggsave(filename = plot_file, plot = plot_trace, width = 9, height = 6)
}

# ---------------------------------------------------
## MODEL COMPARISON (LOO)
# ---------------------------------------------------
loos <- vector("list", length(models))
for (i in seq_along(models)) {
  load(paste0('./results/loo/', exp, '/loo_', models[i], '_', exp, '.rdata'))
  loos[[i]] <- loo
}
names(loos) <- models

loo_comparison <- loo_compare(loos)
print(loo_comparison)

html_content <- kable(loo_comparison, format = "html", table.attr = "class='table table-bordered'")
writeLines(html_content, paste0('./results/summary/', exp, '/loo_report_', exp, '.html'))

df_loo <- tibble(
  model = rownames(loo_comparison),
  elpd_diff = loo_comparison[, "elpd_diff"],
  se_elpd_diff = loo_comparison[, "se_diff"]
)

plot_loo <- ggplot(df_loo, aes(x = reorder(model, elpd_diff), y = elpd_diff, fill = model)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_errorbar(aes(ymin = elpd_diff - se_elpd_diff, ymax = elpd_diff + se_elpd_diff),
                width = 0.2, color = "black") +
  labs(
    title = "Model Comparison (report): ELPD Difference with Error Bars",
    x = "Model",
    y = "ELPD Difference"
  ) +
  theme_minimal()

ggsave(filename = paste0('./results/plots/', exp, '/loo_report_', exp, '.pdf'),
       plot = plot_loo, width = 7, height = 5)

# ---------------------------------------------------
## MODELS SUMMARIES (mu_* + mu_pr[*])
# ---------------------------------------------------
all_summaries <- vector("list", length(models))
for (i in seq_along(models)) {
  fit <- fits[[i]]
  posterior_df <- as_draws_df(fit$draws())

  keep <- grepl('^mu_', colnames(posterior_df)) | grepl('^mu_pr\\[', colnames(posterior_df))
  selected_params <- posterior_df[, keep, drop = FALSE]

  summary_stats <- posterior_summary(selected_params)
  summary_df <- as.data.frame(summary_stats)
  all_summaries[[i]] <- summary_df

  html_file  <- paste0("./results/summary/", exp, "/summary_mu_", models[i], "_", exp, ".html")
  latex_file <- paste0("./results/summary/", exp, "/summary_mu_", models[i], "_", exp, ".tex")

  html_content  <- kable(summary_df, format = "html",  table.attr = "class='table table-bordered'")
  latex_content <- kable(summary_df, format = "latex", booktabs = TRUE)

  writeLines(html_content, html_file)
  writeLines(latex_content, latex_file)
}

################################################################################
##                 PROBA TRANSFORMATION (REPORT MODELS)
################################################################################

invlogit <- function(x) 1 / (1 + exp(-x))
clamp_prob <- function(p) pmin(pmax(p, 1e-6), 1 - 1e-6)

# linear_report: m(p)=alpha*logit(p)+beta  => f(p)=inv_logit(m(p))
f_linear <- function(p, alpha, beta) invlogit(alpha * qlogis(p) + beta)

# BLO_report(+noise): f(p)=inv_logit( w(p)*Lambda(p) + (1-w(p))*anchor )
V <- function(p) p * (1 - p)
w_blo <- function(p, kappa) 1 / (1 + kappa * V(p))

Gamma_clip <- function(p, Delta_minus, Delta_plus) {
  l <- qlogis(clamp_prob(p))
  pmin(pmax(l, Delta_minus), Delta_plus)
}

# FIX: Both BLO_report.stan and BLO_report_noise.stan use
#   denom = sqrt(square(Delta_plus - Delta_minus) + 1e-12)
# The old R code wrongly used (Delta_plus - Delta_minus) without sqrt.
Lambda_map_report <- function(p, Delta_minus, Delta_plus, Psi) {
  G <- Gamma_clip(p, Delta_minus, Delta_plus)
  denom <- sqrt((Delta_plus - Delta_minus)^2 + 1e-12)
  (2 * Psi / denom) * (G - 0.5 * (Delta_plus + Delta_minus))
}

f_blo_pm_report <- function(p, Delta_plus, Delta_minus, Psi, kappa, anchor) {
  pp <- clamp_prob(p)
  wp <- w_blo(pp, kappa)
  L  <- Lambda_map_report(pp, Delta_minus, Delta_plus, Psi)
  invlogit(wp * L + (1 - wp) * anchor)
}

p_grid <- (1:99) / 100

# --- helpers ---
get_est <- function(sum_df, name) {
  if (!name %in% rownames(sum_df)) return(NA_real_)
  sum_df[name, "Estimate"]
}
get_first <- function(sum_df, candidates, default = NA_real_) {
  nm <- candidates[candidates %in% rownames(sum_df)]
  if (length(nm) == 0) return(default)
  sum_df[nm[1], "Estimate"]
}

Phi_approx_R <- function(x) pnorm(x)
safe_exp <- function(x) exp(pmin(x, 20))

# Group params for linear_report:
# Stan exports mu_alpha, mu_beta, mu_lambda, mu_sigma_rep
get_group_linear_report <- function(sum_df) {
  alpha <- get_est(sum_df, "mu_alpha")
  beta  <- get_est(sum_df, "mu_beta")
  if (!is.na(alpha) && !is.na(beta)) return(list(alpha = alpha, beta = beta))

  # fallback from mu_pr[] in linear_report.stan
  # Stan transform_group_means:
  #   out[2] = 2 * Phi_approx(mu_pr[2])   -> alpha
  #   out[3] = mu_pr[3]                    -> beta  (Gaussian, unbounded)
  mu2 <- get_est(sum_df, "mu_pr[2]")
  mu3 <- get_est(sum_df, "mu_pr[3]")
  if (any(is.na(c(mu2, mu3)))) stop("linear_report: cannot find mu_alpha/mu_beta nor mu_pr[2:3].")
  list(
    alpha = 2 * Phi_approx_R(mu2),
    beta  = mu3    # FIX: was -2 + 4*Phi(...); Stan just uses mu_pr[3] directly
  )
}

# Group params for BLO_report models:
# Stan exports mu_Delta_plus, mu_Delta_minus, mu_Psi, mu_kappa, mu_l_anchor
get_group_blo_report <- function(sum_df) {
  Dp   <- get_est(sum_df, "mu_Delta_plus")
  Dm   <- get_est(sum_df, "mu_Delta_minus")
  Psi  <- get_est(sum_df, "mu_Psi")
  kap  <- get_est(sum_df, "mu_kappa")
  anch <- get_est(sum_df, "mu_l_anchor")
  if (!any(is.na(c(Dp, Dm, Psi, kap, anch)))) {
    return(list(Delta_plus = Dp, Delta_minus = Dm, Psi = Psi, kappa = kap, anchor = anch))
  }

  # FIX: fallback rewritten to match actual Stan transform_group_means:
  #   g[2] = exp(min(mu_pr[2], 20))   -> Delta_plus
  #   g[3] = -g[2]                    -> Delta_minus
  #   g[4] = exp(min(mu_pr[4], 20))   -> Psi
  #   g[5] = Phi(mu_pr[5]) * 20       -> kappa
  #   g[6] = mu_pr[6]                 -> l_anchor
  # (mu_pr[3] is a ghost — unused because Delta_minus = -Delta_plus)
  mu2 <- get_est(sum_df, "mu_pr[2]")
  mu4 <- get_est(sum_df, "mu_pr[4]")
  mu5 <- get_est(sum_df, "mu_pr[5]")
  mu6 <- get_est(sum_df, "mu_pr[6]")
  if (any(is.na(c(mu2, mu4, mu5))))
    stop("BLO_report*: cannot find group params (mu_Delta_* or mu_pr[2,4,5,6]).")

  Dp <- safe_exp(mu2)
  list(
    Delta_plus  = Dp,
    Delta_minus = -Dp,
    Psi    = safe_exp(mu4),
    kappa  = Phi_approx_R(mu5) * 20,   # FIX: was safe_exp; Stan uses Phi*20
    anchor = ifelse(is.na(mu6), 0.0, mu6)
  )
}

# --- individual params extraction ---
extract_params_indiv_wide <- function(fit) {
  draws_df <- tryCatch(as_draws_df(fit$draws()), error = function(e) NULL)
  if (is.null(draws_df)) return(NULL)

  pvars <- names(draws_df)[grepl("^params\\[[0-9]+,[0-9]+\\]$", names(draws_df))]
  if (length(pvars) == 0) return(NULL)

  med <- vapply(pvars, function(v) median(draws_df[[v]], na.rm = TRUE), numeric(1))
  idx <- stringr::str_match(pvars, "^params\\[([0-9]+),([0-9]+)\\]$")
  subj <- as.integer(idx[, 2])
  k    <- as.integer(idx[, 3])

  tibble(subj = subj, k = k, median = as.numeric(med)) %>%
    mutate(k = paste0("k", k)) %>%
    pivot_wider(names_from = k, values_from = median) %>%
    arrange(subj)
}

# fallback reconstructions (ONLY if params[n,k] is absent)
reconstruct_linear_report <- function(fit, N) {
  draws_df <- as_draws_df(fit$draws())
  get_med <- function(v) if (v %in% names(draws_df)) median(draws_df[[v]], na.rm = TRUE) else NA_real_

  mu <- vapply(1:4, function(j) get_med(sprintf("mu_pr[%d]", j)), numeric(1))
  sg <- vapply(1:4, function(j) get_med(sprintf("sigma_pr[%d]", j)), numeric(1))
  if (any(is.na(mu)) || any(is.na(sg))) stop("linear_report: missing mu_pr/sigma_pr draws for fallback.")

  out <- vector("list", N)
  for (n in 1:N) {
    raw <- vapply(1:4, function(j) get_med(sprintf("param_raw[%d,%d]", n, j)), numeric(1))
    if (any(is.na(raw))) stop("linear_report: missing param_raw draws for subj ", n)

    # linear_report.stan transform_params:
    #   p[1] = Phi_approx(mu_pr[1] + sigma_pr[1]*raw[1])       -> lambda
    #   p[2] = 2 * Phi_approx(mu_pr[2] + sigma_pr[2]*raw[2])   -> alpha in [0,2]
    #   p[3] = mu_pr[3] + sigma_pr[3]*raw[3]                    -> beta (Gaussian)
    #   p[4] = 2 * Phi_approx(mu_pr[4] + sigma_pr[4]*raw[4])   -> sigma_rep in (0,2)
    out[[n]] <- tibble(
      subj = n,
      lambda    = Phi_approx_R(mu[1] + sg[1] * raw[1]),
      alpha     = 2 * Phi_approx_R(mu[2] + sg[2] * raw[2]),
      beta      = mu[3] + sg[3] * raw[3],   # FIX: was -2 + 4*Phi(...); Stan is Gaussian
      sigma_rep = 2 * Phi_approx_R(mu[4] + sg[4] * raw[4])
    )
  }
  bind_rows(out)
}

reconstruct_blo_report <- function(fit, N) {
  draws_df <- as_draws_df(fit$draws())
  get_med <- function(v) if (v %in% names(draws_df)) median(draws_df[[v]], na.rm = TRUE) else NA_real_

  mu <- vapply(1:7, function(j) get_med(sprintf("mu_pr[%d]", j)), numeric(1))
  sg <- vapply(1:7, function(j) get_med(sprintf("sigma_pr[%d]", j)), numeric(1))
  if (any(is.na(mu)) || any(is.na(sg))) stop("BLO_report*: missing mu_pr/sigma_pr draws for fallback.")

  out <- vector("list", N)
  for (n in 1:N) {
    raw <- vapply(1:7, function(j) get_med(sprintf("param_raw[%d,%d]", n, j)), numeric(1))
    if (any(is.na(raw))) stop("BLO_report*: missing param_raw draws for subj ", n)

    # BLO_report*.stan transform_params:
    #   p[1] = Phi_approx(mu_pr[1] + sigma_pr[1]*raw[1])             -> lambda
    #   p[2] = exp(min(mu_pr[2] + sigma_pr[2]*raw[2], 20))           -> Delta_plus
    #   p[3] = -p[2]                                                  -> Delta_minus (raw[3] UNUSED)
    #   p[4] = exp(min(mu_pr[4] + sigma_pr[4]*raw[4], 20))           -> Psi
    #   p[5] = Phi_approx(mu_pr[5] + sigma_pr[5]*raw[5]) * 20        -> kappa
    #   p[6] = mu_pr[6] + sigma_pr[6]*raw[6]                         -> l_anchor
    #   p[7] = 2 * Phi_approx(mu_pr[7] + sigma_pr[7]*raw[7])         -> sigma_rep
    Dp <- exp(min(mu[2] + sg[2] * raw[2], 20))

    out[[n]] <- tibble(
      subj = n,
      lambda      = Phi_approx_R(mu[1] + sg[1] * raw[1]),
      Delta_plus  = Dp,
      Delta_minus = -Dp,                                                      # FIX: symmetric, not Delta_c ± Delta_h
      Psi         = exp(min(mu[4] + sg[4] * raw[4], 20)),
      kappa       = Phi_approx_R(mu[5] + sg[5] * raw[5]) * 20,               # FIX: was safe_exp; Stan uses Phi*20
      anchor      = mu[6] + sg[6] * raw[6],
      sigma_rep   = 2 * Phi_approx_R(mu[7] + sg[7] * raw[7])
    )
  }
  bind_rows(out)
}

# ---------------------------------------------------
# Extract individual parameters for each model
# ---------------------------------------------------
params_indiv <- vector("list", length(models))
names(params_indiv) <- models

for (i in seq_along(models)) {
  m <- models[i]
  fit <- fits[[i]]

  wide <- extract_params_indiv_wide(fit)

  if (m == "linear_report") {
    # linear_report.stan: params[n,1]=lambda, [2]=alpha, [3]=beta, [4]=sigma_rep
    if (!is.null(wide) && all(c("k1","k2","k3","k4") %in% names(wide))) {
      params_indiv[[m]] <- wide %>%
        transmute(subj,
                  lambda = k1,
                  alpha = k2,
                  beta = k3,
                  sigma_rep = k4)
    } else {
      params_indiv[[m]] <- reconstruct_linear_report(fit, N)
    }
  } else {
    # BLO_report*.stan: params[n,1]=lambda, [2]=Delta_plus, [3]=Delta_minus, [4]=Psi, [5]=kappa, [6]=l_anchor, [7]=sigma_rep
    if (!is.null(wide) && all(c("k1","k2","k3","k4","k5","k6","k7") %in% names(wide))) {
      params_indiv[[m]] <- wide %>%
        transmute(subj,
                  lambda = k1,
                  Delta_plus = k2,
                  Delta_minus = k3,
                  Psi = k4,
                  kappa = k5,
                  anchor = k6,
                  sigma_rep = k7)
    } else {
      params_indiv[[m]] <- reconstruct_blo_report(fit, N)
    }
  }
}

# ---------------------------------------------------
# Build GROUP + INDIV curves
# ---------------------------------------------------
group_curves <- bind_rows(lapply(seq_along(models), function(i) {
  m <- models[i]
  sum_df <- all_summaries[[i]]

  if (m == "linear_report") {
    gp <- get_group_linear_report(sum_df)
    tibble(model = m, level = "group", subj = NA_integer_, p = p_grid,
           fp = f_linear(p_grid, gp$alpha, gp$beta))
  } else {
    gp <- get_group_blo_report(sum_df)
    tibble(model = m, level = "group", subj = NA_integer_, p = p_grid,
           fp = f_blo_pm_report(p_grid, gp$Delta_plus, gp$Delta_minus, gp$Psi, gp$kappa, gp$anchor))
  }
}))

indiv_curves <- bind_rows(lapply(models, function(m) {
  df <- params_indiv[[m]]

  if (m == "linear_report") {
    df %>%
      select(subj, alpha, beta) %>%
      crossing(p = p_grid) %>%
      mutate(model = m, level = "indiv", fp = f_linear(p, alpha, beta)) %>%
      select(model, level, subj, p, fp)
  } else {
    df %>%
      select(subj, Delta_plus, Delta_minus, Psi, kappa, anchor) %>%
      crossing(p = p_grid) %>%
      mutate(model = m, level = "indiv",
             fp = f_blo_pm_report(p, Delta_plus, Delta_minus, Psi, kappa, anchor)) %>%
      select(model, level, subj, p, fp)
  }
}))

curves_all <- bind_rows(group_curves, indiv_curves)

# ---------------------------------------------------
# PLOTS
# ---------------------------------------------------
plot_group <- curves_all %>%
  filter(level == "group") %>%
  ggplot(aes(x = p, y = fp, color = model)) +
  geom_line(linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_minimal() +
  labs(title = paste0("Report: probability transformation — group level (", exp, ")"),
       x = "p", y = "f(p)")

ggsave(
  filename = paste0("./results/plots/", exp, "/report_proba_transform_group_overlay_", exp, ".pdf"),
  plot = plot_group, width = 7, height = 5
)

for (m in models) {
  plt <- curves_all %>%
    filter(level == "indiv", model == m) %>%
    ggplot(aes(x = p, y = fp, group = subj)) +
    geom_line(alpha = 0.15) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal() +
    labs(title = paste0("Report: probability transformation — individuals (", m, ", ", exp, ")"),
         x = "p", y = "f(p)")

  ggsave(
    filename = paste0("./results/plots/", exp, "/report_proba_transform_indiv_spaghetti_", m, "_", exp, ".pdf"),
    plot = plt, width = 7, height = 5
  )
}
