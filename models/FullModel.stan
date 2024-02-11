// Ordered Logistic Multivariate Random Walk model for all items of SCORAD
// - Measurement error by ordered logistic distribution
// - Latent dynamic by multivariate normal random walk
// - Power prior
// - Calibration
// - Treatment effects
// Optional:
// - Population trend (no damping)

functions {
#include /include/functions_OrderedRW.stan
#include /include/get_ragged_bounds.stan
#include /include/get_ts_length.stan
}

data {

  int<lower = 0> N_pt; // Number of patients
  
  // Distributions
  int<lower = 2> M1; // Upper bound of observations
  int<lower = 1> D1; // Number of signs
  int<lower = 2> M2; // Upper bound of observations
  int<lower = 1> D2; // Number of signs
  int<lower = 1, upper = 2> distribution_id[D1 + D2]; // Indicating whether the d-th item use M1 or M2
  
  // Training
  int<lower = 0> N_obs; // Number of non-missing observations
  int<lower = 1, upper = N_pt> k_obs[N_obs]; // Patient index
  int<lower = 1> t_obs[N_obs]; // Time of observation (from 1 to t_max)
  int<lower = 1, upper = D1 + D2> d_obs[N_obs]; // Sign index
  int<lower = 0> y_obs[N_obs]; // Observations
  // Testing
  int<lower = 0> N_test; // Number of predictions to evaluate
  int<lower = 1, upper = N_pt> k_test[N_test]; // Patient index
  int<lower = 1> t_test[N_test]; // Time of prediction
  int<lower = 1, upper = D1 + D2> d_test[N_test]; // Sign index
  int<lower = 0> y_test[N_test]; // True value
  
  // Options
  int<lower = 0, upper = 1> run; // Switch to evaluate the likelihood
  int<lower = 0, upper = 1> independent_items; // Whether to have diagonal correlation matrices or not
  int<lower = 0, upper = 1> trend_known; // Whether to the trend smoothing parameter is known or not
  vector<lower = 0, upper = 1>[D1 + D2] beta_data[trend_known]; // Smoothing parameter
  int<lower = 0> N_agg; // Number of aggregates to compute
  matrix[D1 + D2, N_agg] agg_weights; // Weights of each item in the aggregate
  
#include /include/data_powerprior.stan
#include /include/data_calibration.stan
#include /include/data_dailytreat.stan
  
  // Priors
  vector<lower = 0>[M1 - 1] prior_delta1[D1];
  vector<lower = 0>[M2 - 1] prior_delta2[D2];
  vector[D1 + D2] prior_sigma_meas[2];
  vector[D1 + D2] prior_sigma_lat[2];
  real prior_Omega; 
  real prior_Omega0;
  vector[D1 + D2] prior_mu_y0[2];
  vector<lower = 0>[D1 + D2] prior_sigma_y0[2];
  real prior_ATE[D1 + D2, 2];
  vector<lower = 0>[D1 + D2] prior_beta[2];
  
  // Recommendations
  int<lower = 0> N_rec; // Number of recommendations
  int<lower = 1, upper = N_pt> k_rec[N_rec]; // Patient for which we make the recommendation
  int<lower = 1> t_rec[N_rec]; // Time at which we make the recommendation (prediction at t + 1)
  int<lower = 0> N_actions; // Number of actions to investigate
  matrix[N_actions, D_treat] actions; // Actions
  
}

transformed data {
  // Dealing with two measurements distribution
  int D = D1 + D2; // Total number of items
  int M[D]; // Maximum for all items
  int size_ct = D1 * M1 + D2 * M2; // Size of ragged ct array
  int id_ct[D, 2]; // Index of the first and last values of ct
  int d_sub[D]; // Index of d in sub-arrays corresponding of size D1 or D2
  // Dealing with ragged time-series
  int t_max[N_pt] = get_ts_length(
    append_array(append_array(append_array(append_array(k_obs, k_test), k_cal), k_treat2), k_rec),
    append_array(append_array(append_array(append_array(t_obs, t_test), t_cal), t_treat2), t_rec)
    ); // Length of each time series
  int N = sum(t_max); // Total number of observations
  int id_ts[N_pt, 2] = get_ragged_bounds(t_max); // Index of first and last observation of each patient time-series
  int idx_obs[N_obs]; // index of non-missing observations
  int idx_test[N_test]; // index of predictions
  int idx_rec[N_rec]; // index of recommendations
  int yc_obs[N_obs]; // Categorical y_obs
#include /include/tdata_decl_calibration.stan
#include /include/tdata_decl_dailytreat.stan

  // Dealing with two measurement distributions
  for (d in 1:D) {
    {
      int i1 = 1;
      int i2 = 1;
      if (distribution_id[d] == 1) {
        M[d] = M1;
        d_sub[d] = i1;
        i1 += 1;
      } else {
        M[d] = M2;
        d_sub[d] = i2;
        i2 += 1;
      }
    }
  }
  id_ct = get_ragged_bounds(M);
  
  // Dealing with ragged time-series
  for (i in 1:N_obs) {
    idx_obs[i] = id_ts[k_obs[i], 1] - 1 + t_obs[i];
  }
  for (i in 1:N_test) {
    idx_test[i] = id_ts[k_test[i], 1] - 1 + t_test[i];
  }
  for (i in 1:N_rec) {
    idx_rec[i] = id_ts[k_rec[i], 1] - 1 + t_rec[i];
  }
  
  // Categorical y
  for (i in 1:N_obs) {
    yc_obs[i] = y_obs[i] + 1;
  }

#include /include/tdata_state_dailytreat.stan
#include /include/tdata_state_calibration.stan

}

parameters {
  // Latent dynamic
  vector[D] eta[N]; // Error term, non-centered parametrisation
  cholesky_factor_corr[D] L_param; // Cholesky decomposition of correlation matrix
  vector<lower = 0>[D] sigma_lat; // Vector of standard deviation
  // Measurement distribution
  vector<lower = 0>[D] sigma_meas; // Equivalent standard deviation (not scale) of logistic distribution
  simplex[M1 - 1] delta1[D1]; // Difference between relative cutpoints
  simplex[M2 - 1] delta2[D2]; // Difference between relative cutpoints
  // Initial condition
  cholesky_factor_corr[D] L0_param; // Cholesky decomposition of initial condition correlation matrix
  vector[D] mu_y0; // Population mean y_lat at t0
  vector<lower = 0>[D] sigma_y0; // Population standard deviation of y_lat at t0
  // Treatment effect
  matrix[D, D_treat] ATE; // Average treatment effect (in percentage of score)
  // Trend
  vector<lower = 0, upper = 1>[D] beta_param[1 - trend_known];
  
#include /include/parameters_calibration.stan
#include /include/parameters_dailytreat.stan
  
}

transformed parameters {
  // Latent dynamic and measurement distribution
  matrix[D, D] L;
  matrix[D, D] L0;
  vector[D] s = sigma_meas * sqrt(3) / pi(); // scale of measurement distribution
  vector[M1] ct1[D1]; // Cutpoints in [0, M] space
  vector[M2] ct2[D2]; // Cutpoints in [0, M] space
  vector[size_ct] ct; // ct1 and ct2 concatenated
  vector[size_ct] z_ct; // Cutpoints in affinity space
  vector[D] y_lat[N]; // Latent score
  vector[D] z_lat[N]; // Latent score in affinity space
  vector[D] y0[N_pt]; // Initial latent score
  // Treatment effect
  matrix[D, D_treat] ATE_abs = ATE .* rep_matrix(to_vector(M), D_treat); // Average treatment effect (in units of the score)
  // Trend
  vector[D] trend[N]; // Trend
  vector[D] beta; // Smoothing parameter for the trend
  
#include /include/tparameters_decl_calibration.stan
#include /include/tparameters_decl_dailytreat.stan
  
  if (independent_items == 0) {
    L = L_param;
    L0 = L0_param;
  } else {
    L = diag_matrix(rep_vector(1, D));
    L0 = diag_matrix(rep_vector(1, D));
  }
  
  // Cutpoints
  for (d in 1:D1) {
    ct1[d] = make_ct(delta1[d]);
  }
  for (d in 1:D2) {
    ct2[d] = make_ct(delta2[d]);
  }
  for (d in 1:D) {
    if (distribution_id[d] == 1) {
      ct[id_ct[d, 1]:id_ct[d, 2]] = ct1[d_sub[d]];
    } else {
      ct[id_ct[d, 1]:id_ct[d, 2]] = ct2[d_sub[d]];
    }
    z_ct[id_ct[d, 1]:id_ct[d, 2]] = ct[id_ct[d, 1]:id_ct[d, 2]] / s[d];
  }
  
#include /include/tparameters_state_calibration.stan
#include /include/tparameters_state_dailytreat.stan
  
  // Trend
  if (trend_known) {
    beta = beta_data[1];
  } else {
    beta = beta_param[1];
  }
  
  // Latent dynamic
  for (k in 1:N_pt) {
    y0[k] = mu_y0 + sigma_y0 .* (L0 * eta[id_ts[k, 1]]);
    y_lat[id_ts[k, 1]] = y0[k];
    trend[id_ts[k, 1]] = rep_vector(0, D);
    for (t in (id_ts[k, 1] + 1):id_ts[k, 2]) {
      y_lat[t] = y_lat[t - 1] + trend[t - 1] + ATE_abs * p_treat[t - 1] + sigma_lat .* (L * eta[t]); // Multivariate Random walk
      trend[t] = beta .* (y_lat[t] - y_lat[t - 1]) + (1 - beta) .* trend[t - 1];
    }
  }
  for (i in 1:N) {
    z_lat[i] = y_lat[i] ./ s;
  }
  
}

model {
#include /include/model_dailytreat.stan
#include /include/model_calibration.stan

  for (i in 1:N) {
    eta[i] ~ std_normal();
  }
  // Priors
  // NB: technically the prior should be raised to the power 1-a0 since the power prior is on likelihood and not the posterior
  // But given approximations, weakly informative priors and a0<<1, I don't think it really matters
  L ~ lkj_corr_cholesky(prior_Omega); // LKJ prior for correlation matrix
  L0 ~ lkj_corr_cholesky(prior_Omega0); // LKJ prior for correlation matrix
  for (d in 1:D1) {
    delta1[d] ~ dirichlet(prior_delta1[d]);
  }
  for (d in 1:D2) {
    delta2[d] ~ dirichlet(prior_delta2[d]);
  }
  sigma_meas ./ to_vector(M) ~ lognormal(prior_sigma_meas[1], prior_sigma_meas[2]);
  sigma_lat ./ to_vector(M) ~ lognormal(prior_sigma_lat[1], prior_sigma_lat[2]);
  mu_y0 ./ to_vector(M) ~ normal(prior_mu_y0[1], prior_mu_y0[2]);
  sigma_y0 ./ to_vector(M) ~ normal(prior_sigma_y0[1], prior_sigma_y0[2]);
  // Priors treatment
  for (d in 1:D) {
    ATE[d] ~ normal(prior_ATE[d, 1], prior_ATE[d, 2]); 
  }
  // Prior trend
  if (!trend_known) {
    beta_param[1] ~ beta(prior_beta[1], prior_beta[2]);
  }

  // Power prior
#include /include/model_powerprior.stan
  
  if (run == 1) {
    // Measurement
    for (i in 1:N_obs) {
      yc_obs[i] ~ ordered_logistic(z_lat[idx_obs[i]][d_obs[i]], segment(z_ct, id_ct[d_obs[i], 1], M[d_obs[i]]));
    }
  }
}

generated quantities {
#include /include/gq_decl_dailytreat.stan
#include /include/gq_decl_calibration.stan

  // Additional parameters
  matrix[D, D] Omega = multiply_lower_tri_self_transpose(L); // Correlation matrix
  matrix[D, D] Sigma_lat = quad_form_diag(Omega, sigma_lat); // Covariance matrix
  matrix[D, D] Omega0 = multiply_lower_tri_self_transpose(L0); // Correlation matrix of initial condition
  vector[D] sigma_tot = sqrt(square(sigma_meas) + square(sigma_lat)); // Total noise std for one-step-ahead prediction
  vector[D] sigma_reltot = sigma_tot ./ to_vector(M); // Normalised sigma_tot
  vector[D] rho2 = square(sigma_meas ./ sigma_tot); // Proportion of measurement noise in total noise
  matrix[N_agg, D_treat] ATE_agg = agg_weights' * ATE_abs; // ATE for aggregates
  // Replications of the scores
  matrix[N, D] y_rep; // Replications (of the entire time-series, not just observations)
  matrix[N, N_agg] agg_rep; // Replications of aggregates
  // Predictions
  real lpd[N_test]; // Log predictive density of predictions
  real y_pred[N_test]; // Predictive sample of y_test
  // Recommendations
  matrix[N_rec, D] y_rec[N_actions];
  matrix[N_rec, N_agg] agg_rec[N_actions];
  
#include /include/gq_state_dailytreat.stan
#include /include/gq_state_calibration.stan

  // Replications of the scores
  for (i in 1:N) {
    for (d in 1:D) {
      y_rep[i, d] = ordered_logistic_rng(z_lat[i][d], segment(z_ct, id_ct[d, 1], M[d])) - 1;
    }
  }
  agg_rep = y_rep * agg_weights;

  // Predictions
  for (i in 1:N_test) {
    y_pred[i] = y_rep[idx_test[i], d_test[i]];
    lpd[i] = ordered_logistic_lpmf(y_test[i] + 1 | z_lat[idx_test[i]][d_test[i]], segment(z_ct, id_ct[d_test[i], 1], M[d_test[i]]));
  }

  // Recommendations
  for (a in 1:N_actions) {
    for (i in 1:N_rec) {
      y_rec[a][i] = multi_normal_cholesky_rng(y_lat[idx_rec[i]] + trend[idx_rec[i]] + ATE_abs * actions[a]', diag_matrix(sigma_lat) * L)'; // Linear predictor
      y_rec[a][i] = y_rec[a][i] ./ s';
      for (d in 1:D) {
        y_rec[a][i][d] = ordered_logistic_rng(y_rec[a][i][d], segment(z_ct, id_ct[d, 1], M[d])) - 1; // Measurement
      }
    }
    agg_rec[a] = y_rec[a] * agg_weights;
  }

}
