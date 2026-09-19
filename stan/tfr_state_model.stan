// ==============================================================================
// tfr_state_model.stan
//
// State-level Bayesian demographic model of Total Fertility Rate (TFR)
//
// This model features a hierarchical structure grouping the 36 states by their
// demographic transition timelines. It implements multiple sigma prior regimes
// (Main, Wide, Tight) controlled by the `sigma_prior_regime` data variable
// to support rigorous sensitivity analysis around measurement and process noise.
//
// The global floor (F) prior is centered at log(1.2), reflecting 
// cross-national fertility-transition plateau evidence.
// ==============================================================================

data {
  int<lower=1> N_srs;
  int<lower=1> N_nfhs;
  int<lower=1> S;

  array[N_srs] real y_srs;
  array[N_srs] int<lower=1, upper=S> state_srs;
  array[N_srs] int<lower=1> time_id_srs;

  array[N_nfhs] real y_nfhs;
  array[N_nfhs] int<lower=1, upper=S> state_nfhs;
  array[N_nfhs] int<lower=1> time_id_nfhs;

  int<lower=1> max_T;
  matrix[S, max_T] dt;
  array[S] int<lower=1> T_s;
  vector[S] tfr_init;
  vector[S] tfr_last;
  int<lower=1> n_fore;

  array[S] int<lower=1, upper=5> group_of_state;

  real log_F_mean;
  real log_F_sd;

  int<lower=1, upper=3> sigma_prior_regime;  // 1=MAIN, 2=WIDE, 3=TIGHT
                                              // now ALSO controls timescale_group's prior spread

  real tight_sigma_stat_sd;  // used only when sigma_prior_regime == 3

  int<lower=1, upper=2> model_type;
}

transformed data {
  // Main-model timescale prior: Lognormal(log 50, 0.5^2)
  real ts_mean = log(50);
  real ts_sd;

  if (sigma_prior_regime == 1) {
    ts_sd = 0.5;    // MAIN
  } else if (sigma_prior_regime == 2) {
    ts_sd = 1.0;    // WIDE  -- doubled spread, moves with WIDE sigma regime
  } else {
    ts_sd = 0.25;   // TIGHT -- halved spread, moves with TIGHT sigma regime
  }
}

parameters {
  vector<lower=0>[5] timescale_group;
  real log_F;
  real<lower=0> sigma_stationary;

  vector[S] x1_raw;
  real<lower=0> sigma_init;
  real<lower=0> sigma_srs;
  real<lower=0> sigma_nfhs;

  matrix[S, max_T] z;
}

transformed parameters {
  vector<lower=0>[5] k_group;
  vector<lower=0>[5] sigma_x_group;
  matrix[S, max_T] x;

  for (g in 1:5) {
    k_group[g] = log(10.0) / timescale_group[g];
    sigma_x_group[g] = sigma_stationary * sqrt(2.0 * k_group[g]);
  }

  for (s in 1:S) {
    x[s, 1] = log(tfr_init[s]) + sigma_init * x1_raw[s];
    int g_idx = group_of_state[s];
    real k_val = k_group[g_idx];
    real sx_val = sigma_x_group[g_idx];

    for (t in 2:T_s[s]) {
      real dtt = fmax(dt[s, t], 1e-6);
      real mean_x = log_F + (x[s, t - 1] - log_F) * exp(-k_val * dtt);
      real var_x = (k_val > 1e-8)
                   ? sx_val^2 * (1.0 - exp(-2.0 * k_val * dtt)) / (2.0 * k_val)
                   : sx_val^2 * dtt;
      x[s, t] = mean_x + sqrt(var_x) * z[s, t];
    }
    for (t in (T_s[s] + 1):max_T) {
      x[s, t] = x[s, T_s[s]];
    }
  }
}

model {
  timescale_group ~ lognormal(ts_mean, ts_sd);   // now moves WITH sigma_prior_regime
  log_F ~ normal(log_F_mean, log_F_sd);

  if (sigma_prior_regime == 1) {
    sigma_stationary ~ exponential(2);
    sigma_init       ~ exponential(2);
    sigma_srs        ~ exponential(2);
    sigma_nfhs       ~ exponential(2);
  } else if (sigma_prior_regime == 2) {
    sigma_stationary ~ student_t(3, 0, 2.5) T[0, ];
    sigma_init       ~ student_t(3, 0, 2.5) T[0, ];
    sigma_srs        ~ student_t(3, 0, 2.5) T[0, ];
    sigma_nfhs       ~ student_t(3, 0, 2.5) T[0, ];
  } else {
    sigma_stationary ~ normal(0, tight_sigma_stat_sd) T[0, ];
    sigma_init       ~ normal(0, 0.5);
    sigma_srs        ~ normal(0, 0.5);
    sigma_nfhs       ~ normal(0, 0.5);
  }

  x1_raw       ~ std_normal();
  to_vector(z) ~ std_normal();

  if (model_type == 2) {
    for (n in 1:N_srs) {
      y_srs[n] ~ normal(x[state_srs[n], time_id_srs[n]], sigma_srs);
    }
    for (n in 1:N_nfhs) {
      y_nfhs[n] ~ normal(x[state_nfhs[n], time_id_nfhs[n]], sigma_nfhs);
    }
  }
}

generated quantities {
  array[N_srs] real y_srs_rep;
  array[N_nfhs] real y_nfhs_rep;

  if (model_type == 2) {
    for (n in 1:N_srs) {
      y_srs_rep[n] = normal_rng(x[state_srs[n], time_id_srs[n]], sigma_srs);
    }
    for (n in 1:N_nfhs) {
      y_nfhs_rep[n] = normal_rng(x[state_nfhs[n], time_id_nfhs[n]], sigma_nfhs);
    }
  } else {
    for (n in 1:N_srs) y_srs_rep[n] = 0.0;
    for (n in 1:N_nfhs) y_nfhs_rep[n] = 0.0;
  }
  matrix[S, n_fore] tfr_fore;
  matrix[S, max_T] tfr_hist_fit;
  real F = exp(log_F);
  vector[5] years_to_90pct_closure = timescale_group;

  for (s in 1:S) {
    for (t in 1:max_T) {
      tfr_hist_fit[s, t] = (t <= T_s[s]) ? exp(x[s, t]) : 0.0;
    }
  }

  for (s in 1:S) {
    real x_curr = (model_type == 1) ? log(tfr_last[s]) : x[s, T_s[s]];
    int g_idx = group_of_state[s];
    real k_val = k_group[g_idx];
    real sx_val = sigma_x_group[g_idx];
    for (h in 1:n_fore) {
      real mean_x = log_F + (x_curr - log_F) * exp(-k_val * 1.0);
      real var_x = (k_val > 1e-8)
                   ? sx_val^2 * (1.0 - exp(-2.0 * k_val * 1.0)) / (2.0 * k_val)
                   : sx_val^2 * 1.0;
      x_curr = mean_x + sqrt(var_x) * normal_rng(0, 1);
      tfr_fore[s, h] = exp(x_curr);
    }
  }
}

