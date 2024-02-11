// Data: power prior

// Set a0 to 0 and historical_* to arbitrary values to remove
real<lower = 0, upper = 1> a0; // Discounting factor
real historical_delta1[D1, M1 - 1, 2]; // Historical mean and sd of delta1
real historical_delta2[D2, M2 - 1, 2]; // Historical mean and sd of delta2
real historical_sigma_meas[D1 + D2, 2]; // Historical mean and sd of sigma_meas
real historical_sigma_lat[D1 + D2, 2]; // Historical mean and sd of sigma_meas
real historical_mu_y0[D1 + D2, 2]; // Historical mean and sd of mu_y0
real historical_sigma_y0[D1 + D2, 2]; // Historical mean and sd of sigma_y0