// Transformed data statement: daily treatment usage

// Fill in treat (-1 means missing)
treat2 = rep_array(-1, N, D_treat);
treat = rep_array(-1, N, D_treat);
// When treat2=0
for (i in 1:N_treat2) {
  idx_treat2[i] = id_ts[k_treat2[i], 1] - 1 + t_treat2[i];
  treat2[idx_treat2[i], d_treat2[i]] = treat2_obs[i];
  if (treat2_obs[i] == 0) {
    for (dt in 0:min(2, idx_treat2[i] - id_ts[k_treat2[i], 1])) {
      treat[idx_treat2[i] - dt, d_treat2[i]] = 0;
    }
  }
}
for (d in 1:D_treat)  {
  for (k in 1:N_pt) {
    // Transition 0 to 1
    for (t in (1 + id_ts[k, 1]):id_ts[k, 2]) {
      if (treat2[t - 1, d] == 0 && treat2[t, d] == 1) {
        treat[t, d] = 1;
      }
    }
    // Transition 1 to 0
    for (t in (3 + id_ts[k, 1]):id_ts[k, 2]) {
      if (treat2[t - 1, d] == 1 && treat2[t, d] == 0) {
        treat[t - 3, d] = 1;
      }
    }
  }
}
