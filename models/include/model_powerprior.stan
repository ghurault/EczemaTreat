// Model: power prior

if (a0 > 0) {
  for (d in 1:D) {
    target += a0 * normal_lpdf(sigma_meas[d] | historical_sigma_meas[d, 1], historical_sigma_meas[d, 2]);
    target += a0 * normal_lpdf(sigma_lat[d] | historical_sigma_lat[d, 1], historical_sigma_lat[d, 2]);
    target += a0 * normal_lpdf(mu_y0[d] | historical_mu_y0[d, 1], historical_mu_y0[d, 2]);
    target += a0 * normal_lpdf(sigma_y0[d] | historical_sigma_y0[d, 1], historical_sigma_y0[d, 2]);
  }
  for (d in 1:D1) {
    for (i in 1:(M1 - 1)) {
      target += a0 * normal_lpdf(delta1[d, i] | historical_delta1[d, i, 1], historical_delta1[d, i, 2]);
    }
  }
  for (d in 1:D2) {
    for (i in 1:(M2 - 1)) {
      target += a0 * normal_lpdf(delta2[d, i] | historical_delta2[d, i, 1], historical_delta2[d, i, 2]);
    }
  }
}
