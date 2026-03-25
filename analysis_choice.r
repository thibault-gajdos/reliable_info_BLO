##rm(list = ls(all = TRUE))  ## efface les données
##source('~/thib/projects/tools/R_lib.r')
##setwd('/Users/thibault/thib/projects/reliable_info/reliable_info_identif/')

library(tidyverse)
library(posterior)
library(bayesplot)
library(loo)
library(knitr)

exp <- '8'
models <- c('linear_choice', 'linear_theta_choice', 'BLO_choice', 'BLO_noise_choice')

##load(paste0("./data/data_list_choice_", exp, ".rdata"))
##data_list <- data_list_choice
##N <- data_list$N

# ---------------------------------------------------
## LOAD MODELS FITS  (needed BEFORE diagnosis)
# ---------------------------------------------------
fits <- vector("list", length(models))
for (i in seq_along(models)) {
    model_name = models[[i]]
    file =   paste0("./results/fits/exp", as.character(exp), "/fit_", model_name, "_exp", as.character(exp), ".rdata")
  load(file)
  fits[[i]] <- fit
}

dir.create(paste0("./results/plots/exp", as.character(exp)), recursive = TRUE, showWarnings = FALSE)

dir.create(paste0("./results/summary/exp",as.character(exp)), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------
## MODEL DIAGNOSIS
# ---------------------------------------------------

## Pairs plots
for (i in 1:length(models)) {
  fit <- fits[[i]]
  posterior_samples <- fit$draws()
  posterior_df <- as_draws_df(posterior_samples)
  selected_params <- posterior_df[, grepl("^mu", colnames(posterior_df)) & !grepl("^mu_pr", colnames(posterior_df))]
  plot <- mcmc_pairs(selected_params)
  plot_file <- paste0("./results/plots/exp", as.character(exp), "/pairs_plot_", models[i], "_exp", exp, ".pdf")
  ggsave(plot, file = plot_file)
}

## Traces
for (i in seq_along(models)) {
  fit <- fits[[i]]
  posterior_df <- as_draws_df(fit$draws())

  keep <- grepl('^mu_', colnames(posterior_df)) | grepl('^mu_pr\\[', colnames(posterior_df))
  selected_params <- posterior_df[, keep, drop = FALSE]

  plot_trace <- mcmc_trace(selected_params)

  plot_file <- paste0("./results/plots/exp", as.character(exp), "/trace_plot_", models[i], "_exp", exp, ".pdf")
  ggsave(filename = plot_file, plot = plot_trace, width = 9, height = 6)
}

# ---------------------------------------------------
## MODEL COMPARISON
# ---------------------------------------------------
loos <- vector("list", length(models))
for (i in seq_along(models)) {
  load(paste0('./results/loo/exp', exp, '/loo_', models[i], '_exp', exp, '.rdata'))
  loos[[i]] <- loo
}
names(loos) <- models

loo_comparison <- loo_compare(loos)
print(loo_comparison)

html_content <- kable(loo_comparison, format = "html", table.attr = "class='table table-bordered'")
writeLines(html_content, paste0('./results/summary/exp', exp, '/loo_exp', exp, '.html'))

## Plot ELPD diffs
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
    title = "Model Comparison: ELPD Difference with Error Bars",
    x = "Model",
    y = "ELPD Difference"
  ) +
  theme_minimal()

ggsave(filename = paste0('./results/plots/exp', exp, '/loo_exp', exp, '.pdf'),
       plot = plot_loo, width = 7, height = 5)

# ---------------------------------------------------
## MODELS SUMMARIES (per-model)
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

  html_file  <- paste0("./results/summary/exp", exp, "/summary_mu_", models[i], "_exo", exp, ".html")
  latex_file <- paste0("./results/summary/exp", exp, "/summary_mu_", models[i], "_exp", exp, ".tex")

  html_content  <- kable(summary_df, format = "html",  table.attr = "class='table table-bordered'")
  latex_content <- kable(summary_df, format = "latex", booktabs = TRUE)

  writeLines(html_content, html_file)
  writeLines(latex_content, latex_file)
}

# ---------------------------------------------------
## JOINT SUMMARY TABLES
# ---------------------------------------------------

# Helper: extract named group-level params from a summary df
extract_named_params <- function(sum_df, param_names) {
  out <- list()
  for (p in param_names) {
    if (p %in% rownames(sum_df)) {
      out[[p]] <- sum_df[p, ]
    }
  }
  if (length(out) == 0) return(NULL)
  df <- do.call(rbind, out)
  df$parameter <- rownames(df)
  df
}

# --- Joint table for LINEAR variants ---
linear_models <- c("linear_choice", "linear_theta_choice")
linear_params <- c("mu_lambda", "mu_alpha", "mu_beta", "mu_theta")

linear_joint <- bind_rows(lapply(linear_models, function(m) {
  idx <- which(models == m)
  sum_df <- all_summaries[[idx]]
  df <- extract_named_params(sum_df, linear_params)
  if (is.null(df)) return(NULL)
  df$model <- m
  df
}))

if (nrow(linear_joint) > 0) {
  linear_table <- linear_joint %>%
    select(model, parameter, Estimate, Est.Error, Q2.5, Q97.5) %>%
    arrange(model, parameter)

  html_content <- kable(linear_table, format = "html", digits = 4,
                        table.attr = "class='table table-bordered'",
                        caption = paste0("Linear models — group-level parameters (exp", exp, ")"))
  latex_content <- kable(linear_table, format = "latex", digits = 4,
                         booktabs = TRUE, row.names = FALSE,
                         caption = paste0("Linear models — group-level parameters (exp", exp, ")"))

  writeLines(html_content, paste0("./results/summary/exp", exp, "/joint_summary_linear_exp", exp, ".html"))
  writeLines(latex_content, paste0("./results/summary/exp", exp, "/joint_summary_linear_exp", exp, ".tex"))
}

# --- Joint table for BLO variants ---
blo_models <- c("BLO_choice", "BLO_noise_choice")
blo_params <- c("mu_lambda_recency", "mu_Delta_plus", "mu_Delta_minus",
                "mu_Psi", "mu_kappa", "mu_l_anchor", "mu_theta_choice")

blo_joint <- bind_rows(lapply(blo_models, function(m) {
  idx <- which(models == m)
  sum_df <- all_summaries[[idx]]
  df <- extract_named_params(sum_df, blo_params)
  if (is.null(df)) return(NULL)
  df$model <- m
  df
}))

if (nrow(blo_joint) > 0) {
  blo_table <- blo_joint %>%
    select(model, parameter, Estimate, Est.Error, Q2.5, Q97.5) %>%
    arrange(model, parameter)

  html_content <- kable(blo_table, format = "html", digits = 4,
                        table.attr = "class='table table-bordered'",
                        caption = paste0("BLO models — group-level parameters (exp", exp, ")"))
  latex_content <- kable(blo_table, format = "latex", digits = 4,
                         booktabs = TRUE, row.names = FALSE,
                         caption = paste0("BLO models — group-level parameters (exp", exp, ")"))

  writeLines(html_content, paste0("./results/summary/exp", exp, "/joint_summary_blo_exp", exp, ".html"))
  writeLines(latex_content, paste0("./results/summary/exp", exp, "/joint_summary_blo_exp", exp, ".tex"))
}

################################################################################
##                 PROBA TRANSFORMATION
################################################################################

invlogit <- function(x) 1 / (1 + exp(-x))

# --- BLO building blocks (match Stan) ---
Gamma_clip <- function(p, Delta_minus, Delta_plus) {
  l <- qlogis(p)
  pmin(pmax(l, Delta_minus), Delta_plus)
}

Lambda_map <- function(p, Delta_minus, Delta_plus, Psi) {
  G <- Gamma_clip(p, Delta_minus, Delta_plus)
  denom <- sqrt((Delta_plus - Delta_minus)^2 + 1e-12)
  (2 * Psi / denom) * (G - 0.5 * (Delta_plus + Delta_minus))
}

w_blo <- function(p, kappa) 1 / (1 + kappa * (p * (1 - p)))

f_linear <- function(p, alpha, beta) invlogit(alpha * qlogis(p) + beta)

f_blo <- function(p, Delta_plus, Delta_minus, Psi, kappa, anchor) {
  wp <- w_blo(p, kappa)
  L  <- Lambda_map(p, Delta_minus, Delta_plus, Psi)
  invlogit(wp * L + (1 - wp) * anchor)
}

p_grid <- (1:99) / 100

# --- helpers to extract group parameters robustly ---
get_est <- function(sum_df, name) {
  if (!name %in% rownames(sum_df)) return(NA_real_)
  sum_df[name, "Estimate"]
}

get_first <- function(sum_df, candidates, default = NA_real_) {
  nm <- candidates[candidates %in% rownames(sum_df)]
  if (length(nm) == 0) return(default)
  sum_df[nm[1], "Estimate"]
}

get_group_blo_params <- function(sum_df) {
  # primary path: named generated quantities from Stan
  Dp <- get_est(sum_df, "mu_Delta_plus")
  Dm <- get_est(sum_df, "mu_Delta_minus")

  if (!is.na(Dp) && !is.na(Dm)) {
    Psi    <- get_est(sum_df, "mu_Psi")
    kappa  <- get_est(sum_df, "mu_kappa")
    anchor <- get_first(sum_df, c("mu_l_anchor", "mu_anchor"), default = 0.0)

    if (is.na(Psi))   stop("Cannot find mu_Psi in summary.")
    if (is.na(kappa)) stop("Cannot find mu_kappa in summary.")

    return(list(Delta_plus = Dp, Delta_minus = Dm,
                Psi = Psi, kappa = kappa, anchor = anchor))
  }

  # fallback: reconstruct from mu_pr[]
  mu2 <- get_est(sum_df, "mu_pr[2]")
  mu4 <- get_est(sum_df, "mu_pr[4]")
  mu5 <- get_est(sum_df, "mu_pr[5]")
  mu6 <- get_est(sum_df, "mu_pr[6]")

  if (is.na(mu2) || is.na(mu4) || is.na(mu5))
    stop("Cannot find BLO group params in summary.")

  Dp     <- exp(min(mu2, 20))
  Dm     <- -Dp
  Psi    <- exp(min(mu4, 20))
  kappa  <- pnorm(mu5) * 20
  anchor <- ifelse(is.na(mu6), 0.0, mu6)

  list(Delta_plus = Dp, Delta_minus = Dm,
       Psi = Psi, kappa = kappa, anchor = anchor)
}

# --- individual params extraction ---
is_linear_model <- function(m) m %in% c("linear_choice", "linear_theta_choice")
is_blo_model    <- function(m) m %in% c("BLO_choice", "BLO_noise_choice")

extract_params_indiv_wide <- function(fit, model_name, N = NULL) {
  draws_df <- tryCatch(as_draws_df(fit$draws()), error = function(e) NULL)

  # Case 1: generated quantities includes params[n,k] in draws
  if (!is.null(draws_df)) {
    pvars <- names(draws_df)[grepl("^params\\[[0-9]+,[0-9]+\\]$", names(draws_df))]
    if (length(pvars) > 0) {
      med <- vapply(pvars, function(v) median(draws_df[[v]], na.rm = TRUE), numeric(1))
      idx <- stringr::str_match(pvars, "^params\\[([0-9]+),([0-9]+)\\]$")
      subj <- as.integer(idx[, 2])
      k    <- as.integer(idx[, 3])

      return(
        tibble(subj = subj, k = k, median = as.numeric(med)) %>%
          mutate(k = paste0("k", k)) %>%
          pivot_wider(names_from = k, values_from = median) %>%
          arrange(subj)
      )
    }
  }

  # Case 2 (fallback): reconstruct from fit summary
  if (is.null(N)) stop("extract_params_indiv_wide(): need N for fallback reconstruction.")

  summ <- fit$summary()
  summ <- as.data.frame(summ)
  rownames(summ) <- summ$variable

  get_med <- function(v) {
    if (!v %in% rownames(summ)) return(NA_real_)
    if ("median" %in% names(summ)) return(as.numeric(summ[v, "median"]))
    if ("Mean" %in% names(summ)) return(as.numeric(summ[v, "Mean"]))
    if ("mean" %in% names(summ)) return(as.numeric(summ[v, "mean"]))
    as.numeric(summ[v, 1])
  }

  Phi_approx_R <- function(x) pnorm(x)
  out <- vector("list", N)

  if (model_name == "linear_choice") {
    n_params <- 5
    mu_pr    <- vapply(1:n_params, function(j) get_med(sprintf("mu_pr[%d]", j)), numeric(1))
    sigma_pr <- vapply(1:n_params, function(j) get_med(sprintf("sigma_pr[%d]", j)), numeric(1))
    if (any(is.na(mu_pr)) || any(is.na(sigma_pr)))
      stop("Fallback: cannot find mu_pr/sigma_pr for linear_choice.")
    for (n in 1:N) {
      raw <- vapply(1:n_params, function(j) get_med(sprintf("param_raw[%d,%d]", n, j)), numeric(1))
      if (any(is.na(raw))) stop(paste0("Missing param_raw[", n, ",*]."))
      out[[n]] <- tibble(subj = n,
        k1 = Phi_approx_R(mu_pr[1] + sigma_pr[1] * raw[1]),
        k2 = 10 * Phi_approx_R(mu_pr[2] + sigma_pr[2] * raw[2]),
        k3 = mu_pr[3] + sigma_pr[3] * raw[3],
        k4 = 1.0, k5 = 0.0)
    }

  } else if (model_name == "linear_theta_choice") {
    n_params <- 4
    mu_pr    <- vapply(1:n_params, function(j) get_med(sprintf("mu_pr[%d]", j)), numeric(1))
    sigma_pr <- vapply(1:n_params, function(j) get_med(sprintf("sigma_pr[%d]", j)), numeric(1))
    if (any(is.na(mu_pr)) || any(is.na(sigma_pr)))
      stop("Fallback: cannot find mu_pr/sigma_pr for linear_choice_theta.")
    for (n in 1:N) {
      raw <- vapply(1:n_params, function(j) get_med(sprintf("param_raw[%d,%d]", n, j)), numeric(1))
      if (any(is.na(raw))) stop(paste0("Missing param_raw[", n, ",*]."))
      out[[n]] <- tibble(subj = n,
        k1 = Phi_approx_R(mu_pr[1] + sigma_pr[1] * raw[1]),
        k2 = 5 * Phi_approx_R(mu_pr[2] + sigma_pr[2] * raw[2]),
        k3 = mu_pr[3] + sigma_pr[3] * raw[3],
        k4 = 5 * Phi_approx_R(mu_pr[4] + sigma_pr[4] * raw[4]))
    }

  } else {
    # BLO models: 7 params
    n_params <- 7
    mu_pr    <- vapply(1:n_params, function(j) get_med(sprintf("mu_pr[%d]", j)), numeric(1))
    sigma_pr <- vapply(1:n_params, function(j) get_med(sprintf("sigma_pr[%d]", j)), numeric(1))
    if (any(is.na(mu_pr)) || any(is.na(sigma_pr)))
      stop("Fallback: cannot find mu_pr/sigma_pr for BLO model.")
    for (n in 1:N) {
      raw <- vapply(1:n_params, function(j) get_med(sprintf("param_raw[%d,%d]", n, j)), numeric(1))
      if (any(is.na(raw))) stop(paste0("Missing param_raw[", n, ",*]."))
      Dp <- exp(min(mu_pr[2] + sigma_pr[2] * raw[2], 20))
      out[[n]] <- tibble(subj = n,
        k1 = Phi_approx_R(mu_pr[1] + sigma_pr[1] * raw[1]),
        k2 = Dp, k3 = -Dp,
        k4 = exp(min(mu_pr[4] + sigma_pr[4] * raw[4], 20)),
        k5 = Phi_approx_R(mu_pr[5] + sigma_pr[5] * raw[5]) * 20,
        k6 = mu_pr[6] + sigma_pr[6] * raw[6],
        k7 = Phi_approx_R(mu_pr[7] + sigma_pr[7] * raw[7]) * 10)
    }
  }

  bind_rows(out)
}

# map params[n,k] -> named params per model
map_indiv_linear <- function(df) {
  if (!all(c("k2", "k3") %in% names(df))) stop("linear: expected k2,k3 for alpha,beta.")
  df %>% mutate(alpha = k2, beta = k3)
}

map_indiv_blo <- function(df, sum_df = NULL) {
  if (!all(c("k2", "k3", "k4", "k5") %in% names(df)))
    stop("BLO: expected k2,k3,k4,k5.")
  df$Delta_plus  <- df$k2
  df$Delta_minus <- df$k3
  df$Psi   <- df$k4
  df$kappa <- df$k5
  if ("k6" %in% names(df)) {
    df$anchor <- df$k6
  } else if (!is.null(sum_df)) {
    df$anchor <- get_first(sum_df, c("mu_l_anchor", "mu_anchor", "mu_pr[6]"), default = 0.0)
  } else {
    df$anchor <- 0.0
  }
  df
}

# ---------------------------------------------------
# Extract individual parameters for each model
# ---------------------------------------------------
params_indiv <- lapply(seq_along(models), function(i) {
  m <- models[i]
  sum_df <- all_summaries[[i]]
  df <- extract_params_indiv_wide(fits[[i]], model_name = m, N = N)

  if (is_linear_model(m)) df <- map_indiv_linear(df)
  if (is_blo_model(m))    df <- map_indiv_blo(df, sum_df)

  df$model <- m
  df
})
names(params_indiv) <- models

# ---------------------------------------------------
# Build GROUP curves
# ---------------------------------------------------
group_curves <- bind_rows(lapply(seq_along(models), function(i) {
  m <- models[i]
  sum_df <- all_summaries[[i]]

  if (is_linear_model(m)) {
    alpha <- get_est(sum_df, "mu_alpha")
    beta  <- get_est(sum_df, "mu_beta")
    if (is.na(alpha) || is.na(beta)) stop(paste0(m, ": missing mu_alpha/mu_beta."))
    tibble(model = m, level = "group", subj = NA_integer_,
           p = p_grid, fp = f_linear(p_grid, alpha, beta))
  } else {
    par <- get_group_blo_params(sum_df)
    tibble(model = m, level = "group", subj = NA_integer_,
           p = p_grid,
           fp = f_blo(p_grid, par$Delta_plus, par$Delta_minus,
                      par$Psi, par$kappa, par$anchor))
  }
}))

# ---------------------------------------------------
# Build INDIVIDUAL curves
# ---------------------------------------------------
indiv_curves <- bind_rows(lapply(models, function(m) {
  df <- params_indiv[[m]]
  if (is_linear_model(m)) {
    df %>% select(subj, model, alpha, beta) %>%
      crossing(p = p_grid) %>%
      mutate(level = "indiv", fp = f_linear(p, alpha, beta)) %>%
      select(model, level, subj, p, fp)
  } else {
    df %>% select(subj, model, Delta_plus, Delta_minus, Psi, kappa, anchor) %>%
      crossing(p = p_grid) %>%
      mutate(level = "indiv",
             fp = f_blo(p, Delta_plus, Delta_minus, Psi, kappa, anchor)) %>%
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
  labs(title = paste0("Probability transformation — group level (exp", exp, ")"),
       x = "p", y = "f(p) = inv_logit(m(p))")

ggsave(
  filename = paste0("./results/plots/exp", exp, "/proba_transform_group_overlay_exp", exp, ".pdf"),
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
    labs(title = paste0("Probability transformation — individuals (", m, ", exp", exp, ")"),
         x = "p", y = "f(p)")

  ggsave(
    filename = paste0("./results/plots/exp", exp, "/proba_transform_indiv_spaghetti_", m, "_exp", exp, ".pdf"),
    plot = plt, width = 7, height = 5
  )
}


