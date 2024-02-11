// Parameters declaration: daily treatment usage

vector[D_treat] eta_p01[N_pt];
vector[D_treat] mu_logit_p01;
vector<lower = 0>[D_treat] sigma_logit_p01;
vector[D_treat] eta_p10[N_pt];
vector[D_treat] mu_logit_p10;
vector<lower = 0>[D_treat] sigma_logit_p10;
