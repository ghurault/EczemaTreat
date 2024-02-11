// Transformed parameters statement: calibration

for (d in 1:D) {
    z_cal_ct[id_ct[d, 1]:id_ct[d, 2]] = ct[id_ct[d, 1]:id_ct[d, 2]] / s_cal[d];
}

for (k in 1:N_pt) {
  for (t in id_ts[k, 1]:id_ts[k, 2]) {
    bias[t] = to_vector(include_bias) .* bias0_abs .* exp(-(t - id_ts[k, 1]) ./ tau_bias);
  }
}

