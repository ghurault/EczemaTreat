// Generated quantities statement: declaration

// Replications of the scores as calibrated measurements
for (i in 1:N) {
  for (d in 1:D) {
    y_cal_rep[i, d] = (y_lat[i][d] + bias[i][d]) / s_cal[d]; // Latent score for calibrated measurement
    y_cal_rep[i, d] = ordered_logistic_rng(y_cal_rep[i, d], segment(z_cal_ct, id_ct[d, 1], M[d])) - 1;
  }
}
agg_cal_rep = y_cal_rep * agg_weights;
