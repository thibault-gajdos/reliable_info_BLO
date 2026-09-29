// BLO_choice_noise_cv.stan
// Identical to BLO_choice_noise.stan but with an is_heldout mask in data.
// Held-out trials are excluded from the model likelihood but their
// log_lik is still computed in generated quantities for scoring.

functions {
  // ── (all functions identical to BLO_choice_noise.stan) ──

  vector transform_params(vector mu_pr, vector sigma_pr, row_vector param_raw_n) {
    vector[7] p;
    p[1] = Phi_approx(mu_pr[1] + sigma_pr[1] * param_raw_n[1]);
    p[2] = exp(fmin(mu_pr[2] + sigma_pr[2] * param_raw_n[2], 20));
    p[3] = -p[2];
    p[4] = exp(fmin(mu_pr[4] + sigma_pr[4] * param_raw_n[4], 20));
    p[5] = Phi_approx(mu_pr[5] + sigma_pr[5] * param_raw_n[5]) * 20;
    p[6] = mu_pr[6] + sigma_pr[6] * param_raw_n[6];
    p[7] = Phi_approx(mu_pr[7] + sigma_pr[7] * param_raw_n[7]) * 5;
    return p;
  }

  vector transform_group_means(vector mu_pr) {
    vector[7] g;
    g[1] = Phi_approx(mu_pr[1]);
    g[2] = exp(fmin(mu_pr[2], 20));
    g[3] = -g[2];
    g[4] = exp(fmin(mu_pr[4], 20));
    g[5] = Phi_approx(mu_pr[5]) * 20;
    g[6] = mu_pr[6];
    g[7] = Phi_approx(mu_pr[7]) * 5;
    return g;
  }

  real clamp_prob(real p) { return fmin(fmax(p, 1e-6), 1 - 1e-6); }

  real w(real p, real kappa) { return inv(1 + kappa * p * (1 - p)); }

  real Gamma(real p, real Delta_minus, real Delta_plus) {
    real l = logit(clamp_prob(p));
    return fmin(fmax(l, Delta_minus), Delta_plus);
  }

  real Lambda(real p, real Delta_minus, real Delta_plus, real Psi) {
    real G     = Gamma(p, Delta_minus, Delta_plus);
    real denom = sqrt(square(Delta_plus - Delta_minus) + 1e-12);
    return (2 * Psi / denom) * (G - 0.5 * (Delta_plus + Delta_minus));
  }

  // If Z ~ N(mu, s^2), then E[inv_logit(Z)] ≈ inv_logit( mu / sqrt(1 + (3/pi^2)*s^2 ) )
  // (probit approximation inv_logit(x) ≈ Phi(x / sqrt(pi^2/3)), applied twice; internal noise sd sigma = 1)
  real approx_logit_arg(real mu, real s) {
    real c = 3 / (pi()^2);
    return mu / sqrt(1 + c * square(s));
  }

  vector compute_trial_stats(int S,
                             array[] int color_data,
                             array[] real proba_data,
                             real lambda_recency,
                             real Delta_minus, real Delta_plus,
                             real Psi, real kappa, real l_anchor) {
    vector[2] evidence = rep_vector(0.0, 2);
    real var_w = 0.0;
    vector[2] out;
    if (S > 0) {
      for (s in 1:S) {
        int  c = color_data[s];
        real p = proba_data[s];
        if (c < 1 || c > 2) continue;
        if (p < 0) continue;
        real pp = clamp_prob(p);
        real wp = w(pp, kappa);
        real L  = Lambda(pp, Delta_minus, Delta_plus, Psi);
        real discount = exp(lambda_recency * (s - S));
        real m = wp * L + (1 - wp) * l_anchor;
        evidence[c] += discount * m;
        var_w += square(discount * wp);
      }
    }
    out[1] = evidence[1] - evidence[2];
    out[2] = var_w;
    return out;
  }

  // ── reduce_sum worker: skips held-out trials ──
  real partial_sum(array[] int slice_indices,
                   int start, int end,
                   vector mu_pr, vector sigma_pr,
                   array[] int Tsubj,
                   array[,] int sample,
                   array[,,] int color,
                   array[,,] real proba,
                   array[,] int choice,
                   array[,] int is_heldout,   // CV
                   matrix param_raw) {
    real lp = 0;
    for (i in 1:size(slice_indices)) {
      int n = slice_indices[i];
      vector[7] params = transform_params(mu_pr, sigma_pr, param_raw[n]);
      real lam   = params[1];
      real Dp    = params[2];
      real Dm    = params[3];
      real Psi   = params[4];
      real kap   = params[5];
      real anch  = params[6];
      real theta = params[7];

      for (t in 1:Tsubj[n]) {
        if (is_heldout[n, t] == 1) continue;   // CV

        int S = sample[n, t];
        if (S < 1) continue;
        int ch = choice[n, t];
        if (ch < 1 || ch > 2) continue;

        vector[2] evidence = rep_vector(0.0, 2);
        real var_w = 0.0;
        for (s in 1:S) {
          int  c = color[n, t, s];
          real p = proba[n, t, s];
          if (c < 1 || c > 2) continue;
          if (p < 0) continue;
          real pp = clamp_prob(p);
          real wp = w(pp, kap);
          real L  = Lambda(pp, Dm, Dp, Psi);
          real discount = exp(lam * (s - S));
          real m = wp * L + (1 - wp) * anch;
          evidence[c] += discount * m;
          var_w += square(discount * wp);
        }
        real l_trial     = evidence[1] - evidence[2];
        real sigma_trial = sqrt(var_w + 1e-8);
        real z_mu = theta * l_trial;
        real z_sd = theta * sigma_trial;
        lp += bernoulli_logit_lpmf(ch == 1 | approx_logit_arg(z_mu, z_sd));
      }
    }
    return lp;
  }
}

data {
  int<lower=1> N;
  int<lower=1> T_max;
  int<lower=1> I_max;
  array[N] int<lower=1> Tsubj;
  array[N, T_max] int<lower=-1> sample;
  array[N, T_max, I_max] int<lower=-1, upper=2> color;
  array[N, T_max, I_max] real<lower=-1, upper=1> proba;
  array[N, T_max] int<lower=-1, upper=2> choice;
  int<lower=5> grainsize;

  // ── NEW: 1 = held out from training likelihood ──
  array[N, T_max] int<lower=0, upper=1> is_heldout; // CV
}

parameters {
  vector[7] mu_pr;
  vector<lower=0, upper=10>[7] sigma_pr;
  matrix[N, 7] param_raw;
}

model {
  mu_pr ~ std_normal();
  sigma_pr ~ normal(0, 0.5);
  to_vector(param_raw) ~ std_normal();

  array[N] int indices;
  for (n in 1:N) indices[n] = n;

  target += reduce_sum(partial_sum, indices, grainsize,
                       mu_pr, sigma_pr,
                       Tsubj, sample, color, proba, choice,
                       is_heldout,          // CV
                       param_raw);
}

generated quantities {
  matrix[N, 7] params;
  array[N, T_max] int y_pred = rep_array(-1, N, T_max);
  array[N, T_max] real pred_proba    = rep_array(-1.0, N, T_max);
  array[N, T_max] real diff_evidence = rep_array(-1.0, N, T_max);
  vector[sum(Tsubj)] log_lik;

  vector[7] g = transform_group_means(mu_pr);
  real mu_lambda_recency = g[1];
  real mu_Delta_plus     = g[2];
  real mu_Delta_minus    = g[3];
  real mu_Psi            = g[4];
  real mu_kappa          = g[5];
  real mu_l_anchor       = g[6];
  real mu_theta_choice   = g[7];

  int k = 0;
  for (n in 1:N) {
    params[n] = (transform_params(mu_pr, sigma_pr, param_raw[n]))';

    for (t in 1:Tsubj[n]) {
      k += 1;
      int S  = sample[n, t];
      int ch = choice[n, t];

      if (S < 1 || ch < 1 || ch > 2) {
        log_lik[k] = 0;
        continue;
      }

      array[I_max] int  color_trial;
      array[I_max] real proba_trial;
      for (i in 1:I_max) {
        color_trial[i] = color[n, t, i];
        proba_trial[i] = proba[n, t, i];
      }

      vector[2] stats = compute_trial_stats(
        S, color_trial, proba_trial,
        params[n, 1], params[n, 3], params[n, 2],
        params[n, 4], params[n, 5], params[n, 6]
      );
      real l_trial = stats[1];
      diff_evidence[n, t] = l_trial;

      real logit_arg = approx_logit_arg(
        params[n, 7] * l_trial,
        params[n, 7] * sqrt(stats[2] + 1e-8)
      );

      real p_blue = inv_logit(logit_arg);
      y_pred[n, t]    = bernoulli_rng(p_blue) ? 1 : 2;
      pred_proba[n, t] = p_blue;
      log_lik[k]      = bernoulli_logit_lpmf(ch == 1 | logit_arg);
      // log_lik is computed for ALL trials (train + held-out)
      // Filter in R using the is_heldout mask
    }
  }
}
