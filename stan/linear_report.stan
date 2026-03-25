//*  ---------------------------------------------------
//*               FUNCTIONS
//*  ---------------------------------------------------

functions {

  // ===== Transform parameters =====
  vector transform_params(vector mu_pr, vector sigma_pr, row_vector param_raw_n) {
    vector[4] p;

    // 1 = lambda_recency in (0,1)
    // 2 = alpha in [0,2]
    // 3 = beta in [-2,2]
    // 4 = sigma_rep in (0,2)
    p[1] = Phi_approx(mu_pr[1] + sigma_pr[1] * param_raw_n[1]);               // lambda_recency
    p[2] = 2 * Phi_approx(mu_pr[2] + sigma_pr[2] * param_raw_n[2]);           // alpha
    p[3] = mu_pr[3] + sigma_pr[3] * param_raw_n[3];      // beta
    p[4] = 2 * Phi_approx(mu_pr[4] + sigma_pr[4] * param_raw_n[4]);           // sigma_rep

    return p;
  }

  // ===== Transform group means =====
  vector transform_group_means(vector mu_pr) {
    vector[4] out;

    out[1] = Phi_approx(mu_pr[1]);         // mu_lambda
    out[2] = 2 * Phi_approx(mu_pr[2]);     // mu_alpha
    out[3] = mu_pr[3];// mu_beta
    out[4] = 2 * Phi_approx(mu_pr[4]);     // mu_sigma_rep

    return out;
  }

  // ===== bounded log-odds linear mapping =====
  real m_linear(real p, real alpha, real beta) {
    return alpha * logit(p) + beta;
  }

  // ===== Worker function for reduce_sum =====
  real partial_sum(array[] int slice_indices,
                   int start, int end,
                   vector mu_pr, vector sigma_pr,
                   array[] int Tsubj,
                   array[,] int sample,
                   array[,,] int color,
                   array[,,] real proba,
                   array[,] real lo_proba,
                   matrix param_raw) {
    real lp = 0;

    for (i in 1:size(slice_indices)) {
      int n = slice_indices[i];

      // Subject-level parameters
      vector[4] params = transform_params(mu_pr, sigma_pr, param_raw[n]);

      for (t in 1:Tsubj[n]) {
        vector[2] evidence = rep_vector(0.0, 2);

        int S = sample[n, t];
        if (S < 1) continue;

        for (s in 1:S) {
          int c  = color[n, t, s];
          real p = proba[n, t, s];

          if (c < 1 || c > 2) continue;
          if (p < 0) continue; // padding is -1

          // extra safety
          if (is_nan(p) || is_inf(p)) continue;
          if (p <= 0 || p >= 1) continue;

          real m = m_linear(p, params[2], params[3]);
          real discount = exp(params[1] * (s - S));
          real add = discount * m;

          if (!(is_nan(add) || is_inf(add))) evidence[c] += add;
        }

        if (is_nan(evidence[1]) || is_nan(evidence[2]) || is_inf(evidence[1]) || is_inf(evidence[2]))
          continue;

        real l_trial = evidence[1] - evidence[2];

        // ONLY report noise
        real sigma_trial = sqrt(square(params[4]) + 1e-8);

        lp += normal_lpdf(lo_proba[n, t] | l_trial, sigma_trial);
      }
    }

    return lp;
  }

  vector compute_evidence_lin(int sample_size,
                              array[] int color_data,
                              array[] real proba_data,
                              real lambda_recency,
                              real alpha,
                              real beta) {
    vector[2] evidence = rep_vector(0.0, 2);

    for (s in 1:sample_size) {
      int c  = color_data[s];
      real p = proba_data[s];

      if (c < 1 || c > 2) continue;
      if (p < 0) continue;
      if (is_nan(p) || is_inf(p)) continue;
      if (p <= 0 || p >= 1) continue;

      real m = m_linear(p, alpha, beta);
      real discount = exp(lambda_recency * (s - sample_size));
      real add = discount * m;
      if (!(is_nan(add) || is_inf(add))) evidence[c] += add;
    }

    return evidence;
  }

  real compute_log_lik_lin(int sample_size,
                           array[] int color_data,
                           array[] real proba_data,
                           real lo_proba_obs,
                           real lambda_recency,
                           real alpha,
                           real beta,
                           real sigma_rep) {
    vector[2] evidence = compute_evidence_lin(sample_size, color_data, proba_data,
                                              lambda_recency, alpha, beta);

    if (is_nan(evidence[1]) || is_nan(evidence[2]) || is_inf(evidence[1]) || is_inf(evidence[2]))
      return negative_infinity();

    real l_trial = evidence[1] - evidence[2];
    real sigma_trial = sqrt(square(sigma_rep) + 1e-8);
    return normal_lpdf(lo_proba_obs | l_trial, sigma_trial);
  }
}

//*  ---------------------------------------------------
//*                 DATA
//*  ---------------------------------------------------

data {
  int<lower=1> N;
  int<lower=1> T_max;
  int<lower=1> I_max;

  array[N] int<lower=1> Tsubj;
  array[N, T_max] int<lower=-1> sample;
  array[N, T_max, I_max] int<lower=-1, upper=2> color;
  array[N, T_max, I_max] real<lower=-1, upper=1> proba;

  array[N, T_max] real<lower=-7, upper=7> lo_proba;

  int<lower=5> grainsize;
}

//*  ---------------------------------------------------
//*               PARAMETERS
//*  ---------------------------------------------------

parameters {
  vector[4] mu_pr;
  vector<lower=0>[4] sigma_pr;
  matrix[N, 4] param_raw;
}

//*  ---------------------------------------------------
//*                 MODEL
//*  ---------------------------------------------------

model {
  mu_pr ~ std_normal();
  sigma_pr ~ normal(0, 0.5);
  to_vector(param_raw) ~ std_normal();

  array[N] int indices;
  for (n in 1:N) indices[n] = n;

  target += reduce_sum(partial_sum, indices, grainsize,
                       mu_pr, sigma_pr,
                       Tsubj, sample, color, proba, lo_proba, param_raw);
}

//*  ---------------------------------------------------
//*           GENERATED QUANTITIES
//*  ---------------------------------------------------

generated quantities {
  matrix[N, 4] params;
  array[N, T_max] real y_pred = rep_array(-1.0, N, T_max);
  vector[sum(Tsubj)] log_lik;

  vector[4] mu_tr = transform_group_means(mu_pr);
  real mu_lambda     = mu_tr[1];
  real mu_alpha      = mu_tr[2];
  real mu_beta       = mu_tr[3];
  real mu_sigma_rep  = mu_tr[4];

  int k = 0;
  for (n in 1:N) {
    vector[4] p = transform_params(mu_pr, sigma_pr, param_raw[n]);
    for (j in 1:4) params[n, j] = p[j];

    for (t in 1:Tsubj[n]) {
      k += 1;

      if (sample[n, t] < 1) {
        y_pred[n, t] = -1.0;
        log_lik[k] = 0;
        continue;
      }

      array[I_max] int color_trial;
      array[I_max] real proba_trial;
      for (i in 1:I_max) {
        color_trial[i] = color[n, t, i];
        proba_trial[i] = proba[n, t, i];
      }

      vector[2] evidence = compute_evidence_lin(sample[n, t], color_trial, proba_trial,
                                                p[1], p[2], p[3]);

      if (is_nan(evidence[1]) || is_nan(evidence[2]) || is_inf(evidence[1]) || is_inf(evidence[2])) {
        y_pred[n, t] = -1.0;
        log_lik[k] = negative_infinity();
        continue;
      }

      real mean_trial = evidence[1] - evidence[2];
      real sigma_trial = sqrt(square(p[4]) + 1e-8);

      y_pred[n, t] = normal_rng(mean_trial, sigma_trial);

      log_lik[k] = compute_log_lik_lin(sample[n, t], color_trial, proba_trial,
                                       lo_proba[n, t],
                                       p[1], p[2], p[3],
                                       p[4]);
    }
  }
}
