# Notes -------------------------------------------------------------------

# Fit multivariate model

# Initialisation ----------------------------------------------------------

rm(list = ls()) # Clear Workspace (better to restart the session)

set.seed(2021) # Reproducibility (Stan use a different seed)

source(here::here("analysis", "00_init.R")) # Load libraries, variables and functions

score <- "SCORAD"
dataset <- "PFDC"

#### OPTIONS
model <- ScoradPred(a0 = 0.04, # 0.04
                    independent_items = FALSE,
                    include_calibration = TRUE,
                    include_treatment = TRUE,
                    treatment_names = c("localTreatment", "emollientCream"),
                    include_trend = FALSE,
                    include_recommendations = FALSE)

run <- FALSE
n_chains <- 4
n_it <- 2000
####

stopifnot(
  is_scalar_logical(run),
  is_scalar_wholenumber(n_chains),
  n_chains > 0,
  is_scalar_wholenumber(n_it),
  n_it > 0
)

## Files
file_dict <- get_results_files(outcome = score,
                               model = model$name,
                               dataset = dataset,
                               root_dir = here())

## Parameters
param <- list_parameters(model)
param2 <- list_parameters(model, full_names = TRUE)

# Data --------------------------------------------------------------------

l <- load_PFDC()

# Prepare POSCORAD time-series
POSCORAD <- l$POSCORAD %>%
  rename(Time = Day)
df <- POSCORAD %>%
  select(one_of("Patient", "Time", model$item_spec$Label)) %>%
  pivot_longer(cols = all_of(model$item_spec$Label), names_to = "Label", values_to = "Score") %>%
  drop_na()  %>%
  left_join(model$item_spec[, c("Label", "ItemID")], by = c("Label"))
train <- df %>%
  mutate(Resolution = case_when(Label %in% detail_POSCORAD("Subjective symptoms")$Label ~ 0.1,
                                TRUE ~ 1),
         Score = round(Score / Resolution)) %>%
  select(-Resolution)

# Prepare SCORAD calibration data
if (model$include_calibration) {
  scorad <- l$SCORAD %>%
    rename(Time = Day) %>%
    select(one_of("Patient", "Time", model$item_spec$Label)) %>%
    pivot_longer(cols = all_of(model$item_spec$Label), names_to = "Label", values_to = "Score") %>%
    drop_na()  %>%
    left_join(model$item_spec[, c("Label", "ItemID")], by = c("Label")) 
  cal <- scorad %>%
    mutate(Resolution = case_when(Label %in% detail_POSCORAD("Subjective symptoms")$Label ~ 0.1,
                                  TRUE ~ 1),
           Score = round(Score / Resolution)) %>%
    select(-Resolution)
} else {
  cal <- NULL
}

# Prepare treatment data
treatment_lbl <- paste0(model$treatment_names, "WithinThePast2Days")
if (model$include_treatment) {
  treat <- POSCORAD %>%
    select(all_of(c("Patient", "Time", treatment_lbl))) %>%
    pivot_longer(cols = all_of(treatment_lbl), names_to = "Treatment", values_to = "UsageWithinThePast2Days") %>%
    mutate(Treatment = vapply(Treatment, function(x) {which(x == treatment_lbl)}, numeric(1)) %>% as.numeric()) %>%
    drop_na()
} else {
  treat <- NULL
}

# Prepare recommendation data
if (model$include_recommendations) {
  # Dataset where Time correspond to the time when the action is made, and the scores correspond to Time + 1
  df_rec <- POSCORAD %>%
    group_by(Patient) %>%
    filter(Time == max(Time)) %>%
    ungroup() %>%
    mutate(Time = Time - 1) %>%
    mutate(Recommendation = 1:nrow(.))
} else {
  df_rec <- NULL
}

pt <- unique(df[["Patient"]])

# Stan input
data_stan <- c(prefill_standata_FullModel(model),
               prepare_standata(model, train = train, test = NULL, cal = cal, treat = treat, rec = df_rec))

id <- get_index(bind_rows(train, cal, treat, df_rec))
df <- left_join(df, id, by = c("Patient", "Time"))

# Fitting -----------------------------------------------------------------

if (run) {
  cat("Running model:", model$name, "\n")
  fit <- stan(file = model$stanmodel,
              data = data_stan,
              pars = unlist(param),
              iter = n_it,
              chains = n_chains,
              control = list(adapt_delta = .9),
              init = 0)
  saveRDS(fit, file = file_dict$Fit)
  par <- HuraultMisc::summary_statistics(fit, pars = unlist(param))
  saveRDS(par, file = file_dict$FitPar)
}
