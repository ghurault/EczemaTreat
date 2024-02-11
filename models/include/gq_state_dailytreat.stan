// Generated quantities statement: daily treatment usage

// Replications of the treatments
for (d in 1:D_treat) {
  for (i in 1:N) {
    treat_rep[i, d] = bernoulli_rng(p_treat[i][d]);
  }
  
  for (k in 1:N_pt) {
    for (t in id_ts[k, 1]:(id_ts[k, 1] + 1)) {
      if (treat2[t, d] > -1) {
        treat2_rep[t, d] = treat2[t, d];
      } else {
        treat2_rep[t, d] = bernoulli_rng(1 - (1 - ss1[k][d])^3);
      }
    }
    for (t in (2 + id_ts[k, 1]):id_ts[k, 2]) {
      treat2_rep[t, d] = 1 - (1 - treat_rep[t, d]) * (1 - treat_rep[t - 1, d]) * (1 - treat_rep[t - 2, d]);
    }
  }
}