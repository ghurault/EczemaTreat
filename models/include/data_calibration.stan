// Data: calibration

int<lower = 0> N_cal; // Number of observations for calibration
int<lower = 1, upper = N_pt> k_cal[N_cal]; // Patient index
int<lower = 1, upper = D1 + D2> d_cal[N_cal]; // Item index
int<lower = 1> t_cal[N_cal]; // Time index
int<lower = 0> y_cal[N_cal]; // Calibration value
vector<lower = 0>[D1 + D2] precision_cal; // Ratio of calibration std to measurement std
int<lower = 0, upper = 1> include_bias[D1 + D2]; // Whether to include bias or to set it to 0 (later for subjective symptoms)

// Priors
vector[D1 + D2] prior_bias0[2];
vector[D1 + D2] prior_tau_bias[2];
