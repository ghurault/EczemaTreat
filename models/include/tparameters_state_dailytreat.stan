// Transformed parameters statement: daily treatment usage

for (k in 1:N_pt) {
  p01[k] = inv_logit(mu_logit_p01 + sigma_logit_p01 .* eta_p01[k]);
  p10[k] = inv_logit(mu_logit_p10 + sigma_logit_p10 .* eta_p10[k]);
  ss1[k] = p01[k] ./ (p01[k] + p10[k]);
  p00[k] = 1 - p01[k];
  p11[k] = 1 - p10[k];
}
for (d in 1:D_treat) {
  for (k in 1:N_pt) {
    for (t in id_ts[k, 1]:id_ts[k, 2]) {
      if (treat[t, d] == -1) {
        if (t == id_ts[k, 1]) {
          p_treat[t][d] = ss1[k][d]; // Prior for initial condition
        } else {
          p_treat[t][d] = p11[k][d] * p_treat[t - 1][d] + p01[k][d] * (1 - p_treat[t - 1][d]); // Markov Chain prior
        }
      } else {
        p_treat[t][d] = treat[t, d];
      }
    }
  }
}