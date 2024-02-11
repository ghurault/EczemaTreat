// Transformed parameters declaration: calibration

vector[D] s_cal = s .* precision_cal; // scale of measurement distribution
vector[D] bias0_abs = bias0 .* to_vector(M); // Initial bias in original scale
vector[D] bias[N]; // Calibration bias
vector[size_ct] z_cal_ct; // Cutpoints for calibration in affinity space
