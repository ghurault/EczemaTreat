// Model: daily treatment usage

// Priors
for (k in 1:N_pt) {
  eta_p01[k] ~ std_normal();
  eta_p10[k] ~ std_normal();
}
mu_logit_p01 ~ normal(prior_mu_logit_p01[1], prior_mu_logit_p01[2]);
sigma_logit_p01 ~ normal(prior_sigma_logit_p01[1], prior_sigma_logit_p01[2]);
mu_logit_p10 ~ normal(prior_mu_logit_p10[1], prior_mu_logit_p10[2]);
sigma_logit_p10 ~ normal(prior_sigma_logit_p10[1], prior_sigma_logit_p10[2]);

// Likelihood
if (run == 1) {
  for (d in 1:D_treat) {
    for (k in 1:N_pt) {
      // Likelihood for initial condition
      if (treat[id_ts[k, 1], d] != -1) {
        treat[id_ts[k, 1], d] ~ bernoulli(ss1[k][d]);
      }
      for (t in (1 + id_ts[k, 1]):id_ts[k, 2]) {
        if (treat[t, d] == -1) {
          // If treat is missing and treat2[t] == 1, probability that treatment was used at least once during the past two days
          if (treat2[t, d] == 1 && t > 1 + id_ts[k, 1]) {
            1 ~ bernoulli(1 - (1 - p_treat[t][d]) * (1 - p_treat[t - 1][d]) * (1 - p_treat[t - 2][d]));
          }
        } else {
          // When treatment is observed, Markov Chain likelihood
          treat[t, d] ~ bernoulli(p11[k][d] * p_treat[t - 1][d] + p01[k][d] * (1 - p_treat[t - 1][d]));
        }
      }
    }
  }
}