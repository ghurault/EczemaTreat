// Data: daily treatment usage

int<lower = 1> D_treat; // Number of treatments
int<lower = 0> N_treat2; // Number of "treatment used within the past two days" observations
int<lower = 1, upper = N_pt> k_treat2[N_treat2]; // Patient index corresponding to treat2
int<lower = 1> t_treat2[N_treat2]; // Time index corresponding to treat2
int<lower = 1, upper = D_treat> d_treat2[N_treat2]; // Treatment index
int<lower = 0, upper = 1> treat2_obs[N_treat2]; // Observations of "treatment used within the past two days"

// Priors
real prior_mu_logit_p01[2];
real prior_mu_logit_p10[2];
real prior_sigma_logit_p01[2];
real prior_sigma_logit_p10[2];
