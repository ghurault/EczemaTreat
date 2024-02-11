// Model: calibration

// Priors
bias0 ~ normal(prior_bias0[1], prior_bias0[2]);
tau_bias ~ lognormal(prior_tau_bias[1], prior_tau_bias[2]);

// Likelihood
if (run == 1) {
  for (i in 1:N_cal) {
    yc_cal[i] ~ ordered_logistic((y_lat[idx_cal[i]][d_cal[i]] + bias[idx_cal[i]][d_cal[i]]) / s_cal[d_cal[i]], 
                                 segment(z_cal_ct, id_ct[d_cal[i], 1], M[d_cal[i]]));
  }
}
