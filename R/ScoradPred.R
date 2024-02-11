# Model class -------------------------------------------------------------

#' ScoradPred constructor
#' 
#' Multivariate ordered logistic random walk model with two distributions and optional
#' - power prior
#' - trend
#' - calibration
#'
#' @param independent_items Whether to have diagonal correlations matrices or not
#' @param a0 Forgetting factor for the power prior
#' @param include_trend Whether to include trend in the latent dynamics
#' @param include_calibration Whether we use calibration data
#' @param precision_cal Ratio of calibration std to measurement std
#' @param include_treatment Whether to use treatment data
#' @param treatment_names Vector of treatment names
#' @param include_recommendations Whether to make treatment recommendations (`include_treatment` must be TRUE)
#' @param prior Named list of the model's priors. If `NULL`, uses the default prior for the model (see [default_prior()]).
#' 
#' @return Object of class "ScoradPred" and "EczemaModel"
ScoradPred <- function(independent_items = FALSE,
                       a0 = 0,
                       include_trend = FALSE,
                       include_calibration = FALSE,
                       precision_cal = ifelse(include_calibration, 0.5, 1),
                       include_treatment = FALSE,
                       treatment_names = c(),
                       include_recommendations = FALSE,
                       prior = NULL) {
  
  for (x in list(independent_items,
                 include_trend,
                 include_calibration,
                 include_treatment,
                 include_recommendations)) {
    stopifnot(is_scalar(x),
              is.logical(x))
  }
  
  stopifnot(is_scalar(a0),
            is.numeric(a0),
            between(a0, 0, 1),
            is_scalar(precision_cal),
            is.numeric(precision_cal),
            between(precision_cal, 0, 1),
            include_treatment || !include_recommendations)
  
  if (include_treatment) {
    stopifnot(is.vector(treatment_names, mode = "character"),
              length(treatment_names) > 0)
  } else {
    treatment_names <- c()
  }
  
  M1 <- 100
  M2 <- 3
  item_dict <- ScoradPred_items()
  D1 <- item_dict %>% filter(Distribution == 1) %>% nrow()
  D <- nrow(item_dict)
  D2 <- D - D1
  include_powerprior <- (a0 > 0)
  
  mdl_name <- "ScoradPred"
  suff <- list(
    ifelse(include_powerprior, paste0("h", sprintf("%03d", round(a0 * 100))), ""),
    ifelse(!independent_items, "corr", ""),
    ifelse(include_calibration, "cal", ""),
    ifelse(include_treatment, "treat", ""),
    ifelse(include_trend, "trend", "")
  ) %>%
    setdiff("") %>%
    paste(collapse = "+")
  if (nchar(suff) > 0) {
    mdl_name <- paste0(c(mdl_name, suff), collapse = "+")
  }

  x <- structure(
    list(
      name = mdl_name,
      stanmodel = here("models", paste0("FullModel", ".stan")),
      M1 = M1,
      D1 = D1,
      M2 = M2,
      D2 = D2,
      D = D1 + D2,
      discrete = TRUE,
      item_spec = item_dict,
      independent_items = independent_items,
      a0 = a0,
      include_powerprior = include_powerprior,
      include_trend = include_trend,
      include_calibration = include_calibration,
      precision_cal = precision_cal,
      include_treatment = include_treatment,
      treatment_names = treatment_names,
      include_recommendations = include_recommendations
    ),
    class = c("ScoradPred", "EczemaModel")
  )
  
  if (include_recommendations) {
    x$actions <- get_actions(treatment_names)
  }
  
  x$prior <- c(default_prior(x),
               default_prior_calibration(D),
               default_prior_treatment(D),
               default_prior_trend(D))
  x <- replace_prior(x, prior = prior)
  validate_prior(x)
  
  return(x)
  
}

#' Item characteristics of ScoradPred model
#' 
#' - Distribution 1 corresponds to `M1 = 100` (extent and subjective symptoms)
#' - Distribution 2 corresponds to `M2 = 3` (intensity signs)
#'
#' @return Dataframe from `detail_POSCORAD()` with additional columns
#' - `weight_B`
#' - `weight_C`
#' - `weight_oSCORAD`
#' - `weight_SCORAD`
#' - `ItemID`
#' - `Distribution` (1 or 2)
ScoradPred_items <- function() {
  
  item_dict <- detail_POSCORAD("Items")  %>%
    mutate(M = Maximum / Resolution,
           Component = case_when(Name %in% c("extent") ~ "Extent",
                             Name %in% c("itching", "sleep") ~ "Subjective symptoms",
                             Name %in% detail_POSCORAD("Intensity signs")$Name ~ "Intensity signs"),
           weight_B = as.numeric(Component == "Intensity signs"),
           weight_C = as.numeric(Component == "Subjective symptoms"),
           weight_oSCORAD = case_when(Name == "extent" ~ 0.2,
                                      TRUE ~ 0) + 3.5 * weight_B,
           weight_SCORAD = weight_oSCORAD + weight_C,) %>%
    mutate(across(starts_with("weight_"), ~(.x * Resolution)))
  
  scores1 <- detail_POSCORAD(c("extent", "Subjective symptoms"))$Name
  
  item_dict <- item_dict %>%
    mutate(ItemID = 1:nrow(.),
           Distribution = case_when(Name %in% scores1 ~ 1,
                                    TRUE ~ 2))
  
  return(item_dict)
}

# Default priors ----------------------------------------------------------

#' Default prior for a "multivariate" OrderedRW model
#' 
#' @param max_score Maximum value that the scores can take (same for all items)
#' @param D Number of components
#'
#' @return List
#' @export
default_prior_OrderedMRW <- function(max_score, D) {
  
  dp <- EczemaModel("OrderedRW", max_score = max_score) %>%
    default_prior()

  list(
    delta = matrix(2, nrow = D, ncol = max_score - 1),
    sigma_meas = matrix(dp$sigma_meas, nrow = 2, ncol = D, byrow = FALSE),
    sigma_lat = matrix(dp$sigma_lat, nrow = 2, ncol = D, byrow = FALSE),
    Omega = 10,
    Omega0 = 10,
    mu_y0 = matrix(dp$mu_y0, nrow = 2, ncol = D, byrow = FALSE),
    sigma_y0 = matrix(dp$sigma_y0, nrow = 2, ncol = D, byrow = FALSE)
  )
  
}

#' Default prior for ScoradPred model
#' 
#' The priors are the same as "OrderedMRW" models
#'
#' @param model Object of class "ScoradPred"
#' @param ... Arguments to pass to other methods
#'
#' @return List
default_prior.ScoradPred <- function(model, ...) {
  
  prior <- default_prior_OrderedMRW(max_score = 2, D = model$D1 + model$D2) # max_score does not matter here
  prior1 <- default_prior_OrderedMRW(max_score = model$M1, D = model$D1)
  prior2 <-  default_prior_OrderedMRW(max_score = model$M2, D = model$D2)

  prior$delta <- NULL
  prior <- c(prior,
             list(delta1 = prior1$delta,
                  delta2 = prior2$delta))
  
  return(prior)
  
}

#' Default priors for calibration parameters
#' 
#' - `bias0 ~ normal(x1, x2)` (`bias0` is on a normalised scale)
#' - `tau_bias ~ lognormal(x1, x2)`
#' 
#' NB: there is an identifiability issue when `tau_bias << 1` (bias converges to 0),
#' so we must assume that the order of magnitude of `tau_bias` is more than 1 day.
#'
#' @param D Number of items
#'
#' @return Named list
default_prior_calibration <- function(D) {
  list(
    bias0 = matrix(rep(c(0, 0.1), D), nrow = D, byrow = TRUE) %>% t(),
    tau_bias = matrix(rep(c(1.2 * log(10), 0.6 * log(10)), D), nrow = D, byrow = TRUE) %>% t()
  )
}

#' Default priors for treatment parameters
#' 
#' If there are multiple treatments, the priors are the same regardless of the treatment.
#' 
#' - `mu_logit_p10 ~ normal(x1, x2)`
#' - `mu_logit_p01 ~ normal(x1, x2)`
#' - `sigma_logit_p10 ~ normal(x1, x2)`
#' - `sigma_logit_p01 ~ normal(x1, x2)`
#'
#' @param D Number of items
#'
#' @return Named list
default_prior_treatment <- function(D) {
  list(
    mu_logit_p01 = c(0, 1), 
    mu_logit_p10 = c(0, 1),
    sigma_logit_p01 = c(0, 1.5),
    sigma_logit_p10 = c(0, 1.5),
    ATE = matrix(rep(c(0, 0.1), D), nrow = D, byrow = TRUE)
  )
}

#' Default prior for the trend parameters
#' 
#' - `beta ~ beta(x1, x2)`
#' The element `beta` is 
#'
#' @param D Number of items
#'
#' @return Named list with element `beta` that is a matrix with 2 rows (x1 and x2) and D columns
default_prior_trend <- function(D) {
  list(
    beta = matrix(c(1, 3), nrow = 2, ncol = D, byrow = FALSE)
  )
}

# List parameters ---------------------------------------------------------

#' List available parameters for a "multivariate" OrderedRW model
#'
#' @param max_score Maximum values that the scores can take
#' @param D Number of components
#' @param full_names Whether to return the full names of multi-dimensional parameters
#'
#' @return Named list
list_parameters_OrderedMRW <- function(max_score, D, full_names = FALSE) {
  
  stopifnot(HuraultMisc::is_scalar(full_names),
            is.logical(full_names))
  
  out <- list_parameters("OrderedRW")
  out$Population <- c("Omega", "Omega0", out$Population)
  out$PatientTime <- c(out$PatientTime, "ytot_rep")
  
  if (full_names) {
    
    omg <- expand_grid("i" = 1:D, "j" = 1:D) %>%
      filter(.data$i < .data$j) %>%
      mutate(omg = paste0("Omega[", .data$i, ",", .data$j, "]")) %>%
      pull(omg)
    omg0 <- gsub("Omega", "Omega0", omg)
    
    ct <- expand_2d_parname("ct", D, max_score)
    
    out$Population <- setdiff(out$Population, c("Omega", "Omega0", "ct"))
    out$Population <- c(out$Population, omg, omg0, ct)
    
  }
  
  return(out)
  
}

#' List available parameters for ScoradPred model
#'
#' @param model Object of class "ScoradPred"
#' @param full_names Whether to return the full names of multi-dimensional parameters
#' @param ... Arguments to pass to other methods
#'
#' @return Named list
list_parameters.ScoradPred <- function(model, full_names = FALSE, ...) {
  
  param <- list_parameters_OrderedMRW(max_score = 2,
                                      D = model$D1 + model$D2,
                                      full_names = full_names) # max_score does not matter here
  
  param$Population <- setdiff(param$Population, c("delta", "ct"))
  param$Population <- c(param$Population, "ct1", "ct2", "delta1", "delta2", "sigma_reltot")
  param$Test <- setdiff(param$Test, "cum_err")
  param$PatientTime <- c(setdiff(param$PatientTime, c("ytot_rep")), "agg_rep")

  if (full_names) {
    param$Population <- param$Population[!grepl("^ct\\[\\d+,\\d+\\]$", param$Population)]
    ct1 <- expand_2d_parname("ct1", model$D1, model$M1)
    ct2 <- expand_2d_parname("ct2", model$D2, model$M2)
    param$Population <- c(param$Population, ct1, ct2)
  }
  
  if (model$include_trend) {
    param <- merge_lists(list(param, list_parameters_trend()))
  }
  
  if (model$include_calibration) {
    param <- merge_lists(list(param, list_parameters_calibration()))
  }
  
  if (model$include_treatment) {
    param <- merge_lists(list(param, list_parameters_treatment()))
  }
  
  if (model$include_recommendations) {
    param$Misc <- c("y_rec", "agg_rec")
  }

  return(param)
  
}

list_parameters_calibration <- function() {
  list(
    Population = c("bias0", "tau_bias"), # "bias0_abs"
    PatientTime = c("y_cal_rep", "agg_cal_rep")
  )
}

list_parameters_treatment <- function() {
  list(
    Population = c("mu_logit_p01", "mu_logit_p10", "sigma_logit_p01", "sigma_logit_p10", "ATE", "ATE_agg"),
    Patient = c("p01", "p10", "ss1"),
    PatientTime = c("p_treat", "treat_rep", "treat2_rep")
  )
}

list_parameters_trend <- function() {
  list(
    Population = "beta",
    PatientTime = "trend"
  )
}

# Prefill Stan data -------------------------------------------------------

#' Prefill Stan data
#' 
#' Prefill with item spec, prior, power prior, trend
#'
#' @param model Object of class ScoradPred
#' @param file_historical File posterior summary statistics from historical data
#'
#' @return List
prefill_standata_FullModel <- function(model, file_historical = here("data/par_POSCORAD.rds")) {
  
  # Prefill
  out <- list(
    independent_items = as.numeric(model$independent_items),
    M1 = model$M1,
    M2 = model$M2,
    D1 = model$D1,
    D2 = model$D2,
    distribution_id = model$item_spec$Distribution,
    N_agg = 4,
    agg_weights = as.matrix(model$item_spec[c("weight_B", "weight_C", "weight_oSCORAD", "weight_SCORAD")])
  ) %>%
    add_prior(model$prior)
  
  # Power prior
  par_Derexyl <- readRDS(here("data/par_POSCORAD.rds")) %>%
    filter(Distribution == "Posterior - Derexyl") %>%
    select(-Distribution) %>%
    filter(Variable != "ct")
  par_historical <- model$item_spec %>%
    select(Name, Label, ItemID, Distribution) %>%
    right_join(par_Derexyl, by = c("Name" = "Item"))
  power_prior <- process_powerprior_FullModel(par_historical)
  out <- out %>%
    add_prior(power_prior, prefix = "historical_") %>%
    c(list(a0 = model$a0))
  
  # Trend
  out <- c(out,
           prepare_data_trend(D = model$D, include_trend = model$include_trend))
  
  # Calibration
  if (model$include_calibration) {
    include_bias <- as.numeric(!(model$item_spec$Name %in% c("itching", "sleep"))) # don't calibration subjective symptoms
    precision_cal <- include_bias * model$precision_cal + (1 - include_bias) * 1
  } else {
    include_bias <- rep(0, model$D)
    precision_cal <- rep(1, model$D)
  }
  out <- c(out,
           list(include_bias = include_bias,
                precision_cal = precision_cal))
  
  # NB: for treatment, D_treat inputted in prepare_standata
  
  # Recommendations
  if (model$include_recommendations) {
    data_rec <- list(N_actions = nrow(model$actions),
                     actions = as.matrix(model$actions[, model$treatment_names]))
  } else {
    data_rec <- list(N_actions = 0,
                     actions = matrix(numeric(0),
                                      nrow = 0,
                                      ncol = ifelse(model$include_treatment, length(model$treatment_names), 1)))
    # D_treat=1 when include_treatment=FALSE (implemented as treatment never used)
  }
  out <- c(out, data_rec)
  
  return(out)
  
}

#' List of trend inputs to pass to the Stan sampler
#' NB: need to rename as this goes in the prefill
#'
#' @param D Number of items
#' @param include_trend Whether to include trend or not
#'
#' @return List with elements `trend_known` and `beta_data`
prepare_data_trend <- function(D, include_trend = FALSE) {
  
  if (include_trend) {
    out <- list(
      trend_known = 0,
      beta_data = matrix(numeric(0), nrow = 0, ncol = D)
    )
  } else {
    out <- list(
      trend_known = 1,
      beta_data = matrix(rep(0, D), nrow = 1, ncol = D)
    )
  }
  
  return(out)
}

# Prepare Stan data -------------------------------------------------------

#' Prepare the data list to pass to the Stan sampler
#' 
#' This function only prepares the time-series data.
#'
#' @param model Object of class ScoradPred
#' @param train Training dataframe
#' @param test Testing dataframe
#' @param cal Calibration dataframe
#' @param treat Treatment dataframe
#' @param rec Recommendation dataframe
#' @param ... Arguments to pass to other methods
#'
#' @return List to serve as input to the Stan sampler
prepare_standata.ScoradPred <- function(model, train, test = NULL, cal = NULL, treat = NULL, rec = NULL, ...) {
  
  stopifnot(all(train[["ItemID"]] %in% 1:model$D))
  
  tmp_train <- train %>%
    left_join(model$item_spec, by = "ItemID")
  sc1 <- tmp_train %>%
    filter(Distribution == 1) %>%
    pull(Score)
  stopifnot(all(sc1 %in% 0:model$M1))
  sc2 <- tmp_train %>%
    filter(Distribution == 2) %>%
    pull(Score)
  stopifnot(all(sc2 %in% 0:model$M2))

  out <- prepare_data_lgtd(train = train, test = test, max_score = 1 + max(model$M1, model$M2), discrete = TRUE)
  out$M <- NULL
  
  out$d_obs <- train[["ItemID"]]
  out$d_test <- vector()
  
  out$run <- 1 # do that outside this function for consistency?
  
  if (!is.null(test)) {
    stopifnot(all(test[["ItemID"]] %in% 1:model$D))
    out$d_test <- array(test[["ItemID"]])
  }
  
  # Calibration
  if (model$include_calibration) {
    if (is.null(cal)) {
      warning("include_calibration=TRUE but cal is not supplied")
    }
  } else {
    if (!is.null(cal)) {
      warning("cal is supplied but include_calibration=FALSE will take precedent")
    }
    cal <- NULL
  }
  data_cal <- prepare_data_calibration(df = cal)
  if (length(data_cal$t_cal) > 0 && length(out$t_test) > 0 && max(data_cal$t_cal) >= min(out$t_test)) {
    warning("cal may overlap with test")
  }
  out <- c(out, data_cal)
  
  # Recommendations
  if (model$include_recommendations) {
    if (is.null(rec)) {
      warning("include_recommendations=TRUE but rec is not supplied")
    }
  } else {
    if (!is.null(rec)) {
      warning("rec is supplied but include_recommendations=FALSE will take precedent")
    }
    rec <- NULL
  }
  data_rec <- prepare_data_recommendations(rec)
  out <- c(out, data_rec)
  
  # Treatment
  if (model$include_treatment && !is.null(treat)) {
    stopifnot(max(treat[["Treatment"]]) == length(model$treatment_names))
    data_treat <- prepare_data_treatment(treat)
  } else {
    
    if (model$include_treatment && is.null(treat)) {
      warning("include_treatment=TRUE but treat is not supplied")
    }
    if (!model$include_treatment && !is.null(treat)) {
      warning("treat is supplied but include_treatment=FALSE will take precedent")
    }
    
    id <- bind_rows(train, test, cal, treat, rec) %>%
      get_index()
    data_treat <- list(
      N_treat2 = nrow(id),
      D_treat = 1,
      k_treat2 = id[["Patient"]],
      t_treat2 = id[["Time"]], 
      d_treat2 = rep(1, nrow(id)),
      treat2_obs = rep(0, nrow(id))
    )

  }
  out <- c(out, data_treat)

  return(out)
}

#' List of calibration inputs to pass to Stan sampler
#'
#' @param df (Optional) Dataframe with columns Patient, ItemID, Time and Score
#'
#' @return List
prepare_data_calibration <- function(df = NULL) {
  
  if (is.null(df)) {
    
    out <- list(
      N_cal = 0,
      k_cal = vector(),
      d_cal = vector(),
      t_cal = vector(),
      y_cal = vector()
    )
    
  } else {
    
    stopifnot(
      is.data.frame(df),
      all(c("Patient", "ItemID", "Time", "Score") %in% colnames(df)),
      all(vapply(c("Patient", "ItemID", "Time", "Score"), function(x) {is.numeric(df[[x]])}, logical(1))),
      all(is_wholenumber(df[["Patient"]])),
      all(is_wholenumber(df[["ItemID"]])),
      all(is_wholenumber(df[["Time"]]))
    )
    
    out <- list(
      N_cal = nrow(df),
      k_cal = df[["Patient"]],
      d_cal = df[["ItemID"]],
      t_cal = df[["Time"]],
      y_cal = df[["Score"]]
    )
  }
  
  return(out)
  
}

#' List of treatment inputs to pass to the Stan sampler
#'
#' @param df Dataframe with columns Patient, Time, Treatment and UsageWithinThePast2Days
#'
#' @return List
prepare_data_treatment <- function(df) {
  
  cnames <- c("Patient", "Time", "Treatment", "UsageWithinThePast2Days")
  stopifnot(
    is.data.frame(df),
    all(cnames %in% colnames(df)),
    all(vapply(cnames, function(x) {is.numeric(df[[x]])}, logical(1))),
    all(vapply(cnames, function(x) {all(is_wholenumber(df[[x]]))}, logical(1))),
    all(df[["UsageWithinThePast2Days"]] %in% c(0, 1))
  )
  
  list(
    N_treat2 = nrow(df),
    D_treat = max(df[["Treatment"]]),
    k_treat2 = df[["Patient"]],
    t_treat2 = df[["Time"]],
    d_treat2 = df[["Treatment"]],
    treat2_obs = df[["UsageWithinThePast2Days"]]
  )
}

#' List of recommendation inputs to pass to the Stan sampler
#'
#' @param df Dataframe containing (Patient, Time) of recommendations
#'
#' @return List
prepare_data_recommendations <- function(df = NULL) {
  
  if (is.null(df)) {
    out <- list(
      N_rec = 0,
      k_rec = vector(),
      t_rec = vector()
    )
  } else {
    stopifnot(is.data.frame(df),
              all(c("Patient", "Time") %in% colnames(df)),
              all(vapply(c("Patient", "Time"), function(x) {is.numeric(df[[x]])}, logical(1))),
              all(vapply(c("Patient", "Time"), function(x) {all(is_wholenumber(df[[x]]))}, logical(1))))
    
    out <- list(
      N_rec = nrow(df),
      k_rec = df[["Patient"]],
      t_rec = df[["Time"]]
    )
  }
  
  return(out)
}
