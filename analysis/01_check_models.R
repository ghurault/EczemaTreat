# Notes -------------------------------------------------------------------

# Prior predictive check
# NB: assume no power prior, no trend, no calibration data, no treatment data but correlations between signs

# Initialisation ----------------------------------------------------------

rm(list = ls()) # Clear Workspace (better to restart the session)

set.seed(2021) # Reproducibility (Stan use a different seed)

source(here::here("analysis", "00_init.R")) # Load libraries, variables and functions

score <- "SCORAD"

#### OPTIONS
model <- ScoradPred(independent_items = FALSE)

n_pt <- 16
n_dur <- rpois(n_pt, 50)

run_prior <- TRUE
n_chains <- 4
n_it <- 2000
####

stopifnot(
  is_scalar_wholenumber(n_pt),
  n_pt > 0,
  all(is_wholenumber(n_dur)),
  all(n_dur > 0),
  is_scalar_logical(run_prior),
  is_scalar_wholenumber(n_chains),
  n_chains > 0,
  is_scalar_wholenumber(n_it),
  n_it > 0
)

## Files
file_dict <- get_results_files(outcome = score,
                               model = model$name)

if (run_prior) {
  compiled_model <- stan_model(model$stanmodel)
}

## Parameters
param <- list_parameters(model)
param2 <- list_parameters(model, full_names = TRUE)
param[c("PatientTime", "Test")] <- NULL

id <- get_index2(n_dur)

# Prepare Stan input ------------------------------------------------------

l <- make_empty_data(N_patient = n_pt, t_max = n_dur, max_score = max(model$M1, model$M2), discrete = TRUE)
l$Training$ItemID <- 1
l$Testing$ItemID <- 1
data_stan <- prepare_standata(model, train = l$Training, test = l$Testing)
data_stan[c("N_obs", "d_obs", "k_obs", "t_obs", "y_obs", "run")] <- NULL
data_stan <- c(data_stan,
               list(N_obs = 0,
                    d_obs = vector(),
                    k_obs = vector(),
                    t_obs = vector(),
                    y_obs = vector(),
                    run = 0))

data_prior <- c(prefill_standata_FullModel(model),
                data_stan)

# Prior predictive check -------------------------------------------------

if (run_prior) {
  fit_prior <- sampling(compiled_model,
                        data = data_prior,
                        pars = unlist(param),
                        iter = n_it,
                        chains = n_chains)
  saveRDS(fit_prior, file = here(file_dict$PriorFit))
  par0 <- HuraultMisc::summary_statistics(fit_prior, pars = unlist(param))
  saveRDS(par0, file = here(file_dict$PriorPar))
}
