//*  ---------------------------------------------------
//*    BLO report model (NO internal noise)
//*    Symmetric Delta (Delta+ = -Delta-), structured
//*    like BLO_choice_noise.stan
//*  ---------------------------------------------------

//*  ---------------------------------------------------
//*               FUNCTIONS
//*  ---------------------------------------------------

functions {

  // ===== Transform parameters =====
  // 1 = lambda_recency in [0,1]
  // 2 = Delta_plus  (>0, and Delta_minus = -Delta_plus)
  // 3 = Delta_minus (= -Delta_plus, not a free parameter)
  // 4 = Psi > 0
  // 5 = kappa in [0,20]
  // 6 = l_anchor (real)
  // 7 = sigma_rep in (0,2)
  vector transform_params(vector mu_pr, vector sigma_pr, row_vector param_raw_n) {
    vector[7] p;
    p[1] = Phi_approx(mu_pr[1] + sigma_pr[1] * param_raw_n[1]);              // lambda_recency in [0,1]
    p[2] = exp(fmin(mu_pr[2] + sigma_pr[2] * param_raw_n[2], 20));           // Delta_plus > 0
    p[3] = -p[2];                                                             // Delta_minus = -Delta_plus
    p[4] = exp(fmin(mu_pr[4] + sigma_pr[4] * param_raw_n[4], 20));           // Psi > 0
    p[5] = Phi_approx(mu_pr[5] + sigma_pr[5] * param_raw_n[5]) * 20;        // kappa in [0,20]
    p[6] = mu_pr[6] + sigma_pr[6] * param_raw_n[6];                          // l_anchor (gaussian)
    p[7] = 2 * Phi_approx(mu_pr[7] + sigma_pr[7] * param_raw_n[7]);          // sigma_rep in (0,2)
    return p;
  }

  vector transform_group_means(vector mu_pr) {
    vector[7] g;
    g[1] = Phi_approx(mu_pr[1]);              // mu_lambda_recency
    g[2] = exp(fmin(mu_pr[2], 20));           // mu_Delta_plus
    g[3] = -g[2];                              // mu_Delta_minus
    g[4] = exp(fmin(mu_pr[4], 20));           // mu_Psi
    g[5] = Phi_approx(mu_pr[5]) * 20;        // mu_kappa
    g[6] = mu_pr[6];                           // mu_l_anchor
    g[7] = 2 * Phi_approx(mu_pr[7]);          // mu_sigma_rep
    return g;
  }

  // ===== clamping =====
  real clamp_prob(real p) {
    return fmin(fmax(p, 1e-6), 1 - 1e-6);
  }

  // ===== BLO components =====
  real w(real p, real kappa) {
    return inv(1 + kappa * p * (1 - p));
  }

  real Gamma(real p, real Delta_minus, real Delta_plus) {
    real l = logit(clamp_prob(p));
    return fmin(fmax(l, Delta_minus), Delta_plus);
  }

  real Lambda(real p, real Delta_minus, real Delta_plus, real Psi) {
    real G = Gamma(p, Delta_minus, Delta_plus);
    real denom = sqrt(square(Delta_plus - Delta_minus) + 1e-12);
    return (2 * Psi / denom) * (G - 0.5 * (Delta_plus + Delta_minus));
  }

  // ===== Evidence helper =====
  vector compute_evidence_blo(int S,
                              array[] int color_data,
                              array[] real proba_data,
                              real lambda_recency,
                              real Delta_minus, real Delta_plus,
                              real Psi, real kappa, real l_anchor) {
    vector[2] evidence = rep_vector(0.0, 2);

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
      }
    }

    return evidence;
  }

  // ===== reduce_sum worker (report model, NO internal noise) =====
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

      vector[7] params = transform_params(mu_pr, sigma_pr, param_raw[n]);

      // locals for speed
      real lam   = params[1];
      real Dp    = params[2];
      real Dm    = params[3];
      real Psi   = params[4];
      real kap   = params[5];
      real anch  = params[6];
      real sig_r = params[7];

      real sigma_trial = sqrt(square(sig_r) + 1e-8);

      for (t in 1:Tsubj[n]) {
        int S = sample[n, t];
        if (S < 1) continue;

        vector[2] evidence = rep_vector(0.0, 2);

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
        }

        real mean_trial = evidence[1] - evidence[2];

        lp += normal_lpdf(lo_proba[n, t] | mean_trial, sigma_trial);
      }
    }
    return lp;
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
//*              PARAMETERS
//*  ---------------------------------------------------

parameters {
  vector[7] mu_pr;
  vector<lower=0, upper=10>[7] sigma_pr;
  matrix[N, 7] param_raw;
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
//*          GENERATED QUANTITIES
//*  ---------------------------------------------------

generated quantities {
  matrix[N, 7] params;
  array[N, T_max] real y_pred = rep_array(-1.0, N, T_max);
  array[N, T_max] real diff_evidence = rep_array(-1.0, N, T_max);
  vector[sum(Tsubj)] log_lik;

  // --- group-level (interpretable) means ---
  vector[7] g = transform_group_means(mu_pr);
  real mu_lambda_recency = g[1];
  real mu_Delta_plus     = g[2];
  real mu_Delta_minus    = g[3];
  real mu_Psi            = g[4];
  real mu_kappa          = g[5];
  real mu_l_anchor       = g[6];
  real mu_sigma_rep      = g[7];

  int k = 0;
  for (n in 1:N) {
    // subject-level transformed params
    params[n] = (transform_params(mu_pr, sigma_pr, param_raw[n]))';

    for (t in 1:Tsubj[n]) {
      k += 1;

      if (sample[n, t] < 1) {
        y_pred[n, t] = -1.0;
        diff_evidence[n, t] = -1.0;
        log_lik[k] = 0;
        continue;
      }

      // pack trial arrays for helper
      array[I_max] int color_trial;
      array[I_max] real proba_trial;
      for (i in 1:I_max) {
        color_trial[i] = color[n, t, i];
        proba_trial[i] = proba[n, t, i];
      }

      vector[2] evidence = compute_evidence_blo(
        sample[n, t],
        color_trial, proba_trial,
        params[n, 1],                   // lambda_recency
        params[n, 3], params[n, 2],     // Delta_minus, Delta_plus
        params[n, 4],                   // Psi
        params[n, 5],                   // kappa
        params[n, 6]                    // l_anchor
      );

      real mean_trial  = evidence[1] - evidence[2];
      real sigma_trial = sqrt(square(params[n, 7]) + 1e-8);

      diff_evidence[n, t] = mean_trial;
      y_pred[n, t] = normal_rng(mean_trial, sigma_trial);
      log_lik[k] = normal_lpdf(lo_proba[n, t] | mean_trial, sigma_trial);
    }
  }
}
