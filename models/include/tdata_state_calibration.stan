// Transformed data statement: calibration

for (i in 1:N_cal) {
  idx_cal[i] = id_ts[k_cal[i], 1] - 1 + t_cal[i];
}

for (i in 1:N_cal) {
  yc_cal[i] = y_cal[i] + 1;
}
