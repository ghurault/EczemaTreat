# Notes -------------------------------------------------------------------

# Run validation for multivariate models
# NB: use t_horizon=4 for model comparison but t_horizon=1 for recommendations

# Initialisation ----------------------------------------------------------

rm(list = ls()) # Clear Workspace (better to restart the session)

set.seed(2021) # Reproducibility (Stan use a different seed)

source(here::here("analysis", "00_init.R")) # Load libraries, variables and functions
library(foreach)
library(doParallel)

score <- "SCORAD"
dataset <- "PFDC"

#### OPTIONS
model <- ScoradPred(a0 = 0.04, # 0.04
                    independent_items = FALSE,
                    include_calibration = TRUE,
                    include_treatment = TRUE,
                    treatment_names = c("localTreatment", "emollientCream"),
                    include_trend = FALSE,
                    include_recommendations = TRUE)
# set include_recommendations the same as include_treatment

run <- FALSE
t_horizon <- 4
n_chains <- 4
n_it <- 2000
n_cluster <- 4
####

stopifnot(
  is_scalar_logical(run),
  is_scalar_wholenumber(n_chains),
  n_chains > 0,
  is_scalar_wholenumber(n_it),
  n_it > 0,
  is_scalar_wholenumber(t_horizon),
  t_horizon > 0,
  is_scalar_wholenumber(n_cluster),
  between(n_cluster, 1, floor((parallel::detectCores() - 2) / n_chains))
)

## Parameters
param <- c("lpd", "agg_rep", "y_pred")
if (model$include_recommendations) {
  param <- c(param, "y_rec", "agg_rec", "p_treat")
}

## Files
outcomes <- detail_POSCORAD()$Name
# Validation files
file_dict <- lapply(outcomes,
                    function(x) {
                      get_results_files(outcome = x,
                                        model = model$name,
                                        dataset = dataset,
                                        val_horizon = t_horizon,
                                        root_dir = here())
                    })
names(file_dict) <- outcomes
# Recommendation files
rec_files <- get_recommendation_files(outcome = score,
                                      model = model$name,
                                      dataset = dataset,
                                      val_horizon = t_horizon,
                                      root_dir = here())

if (run) {
  compiled_model <- rstan::stan_model(model$stanmodel)
}

# Prepare Stan input ------------------------------------------------------

l <- load_PFDC()

POSCORAD <- l$POSCORAD %>%
  rename(Time = Day)

# Prefill Stan input
data_stan0 <- prefill_standata_FullModel(model)

# Prepare dataset (model$item_spec controls the indexing of items)
df <- POSCORAD %>%
  select(one_of("Patient", "Time", model$item_spec$Label)) %>%
  pivot_longer(cols = all_of(model$item_spec$Label), names_to = "Label", values_to = "Score") %>%
  drop_na() %>%
  left_join(model$item_spec[, c("Label", "ItemID")], by = c("Label")) %>%
  select(-Label)

if (model$include_calibration) {
  # Format SCORAD
  scorad <- l$SCORAD %>%
    rename(Time = Day) %>%
    select(one_of("Patient", "Time", model$item_spec$Label)) %>%
    pivot_longer(cols = all_of(model$item_spec$Label), names_to = "Label", values_to = "Score") %>%
    drop_na()  %>%
    left_join(model$item_spec[, c("Label", "ItemID")], by = c("Label"))
  scorad <- scorad %>% mutate(Iteration = get_fc_iteration(Time, horizon = t_horizon))
}
df <- df %>% mutate(Iteration = get_fc_iteration(Time, horizon = t_horizon))

if (model$include_treatment) {
  
  treatment_lbl <- paste0(model$treatment_names, "WithinThePast2Days")

  treat <- POSCORAD %>%
    select(all_of(c("Patient", "Time", treatment_lbl))) %>%
    pivot_longer(cols = all_of(treatment_lbl), names_to = "Treatment", values_to = "UsageWithinThePast2Days") %>%
    mutate(Treatment = vapply(Treatment, function(x) {which(x == treatment_lbl)}, numeric(1)) %>% as.numeric()) %>%
    drop_na()
  
  treat <- treat %>% mutate(Iteration = get_fc_iteration(Time, horizon = t_horizon))
  
}

# Nothing to prepare for recommendation (or trend)

pt <- unique(df[["Patient"]])
t_max <- df %>%
  group_by(Patient) %>%
  summarise(LastTime = max(Time)) %>%
  ungroup()

# Forward chaining --------------------------------------------------------

train_it <- get_fc_training_iteration(df[["Iteration"]])
fc_it <- detail_fc_training(df, horizon = t_horizon)

if (run) {
  
  cl <- makeCluster(n_cluster, outfile =  "")
  registerDoParallel(cl)
  
  for (j in 1:length(file_dict)) {
    dir.create(file_dict[[j]]$ValDir)
  }
  if (model$include_recommendations) {
    dir.create(rec_files$RecDir)
  }
  
  out <- foreach(i = rev(seq_along(train_it))) %dopar% {
    it <- train_it[i]
    
    # Need to reload functions and libraries
    source(here::here("analysis", "00_init.R"))
    
    duration <- Sys.time()
    cat(glue::glue("Starting iteration {it}"), sep = "\n")
    
    ####
    
    # Split dataset
    split <- lapply(1:model$D,
                    function(d) {
                      df %>%
                        filter(ItemID == d) %>%
                        select(-ItemID) %>%
                        split_fc_dataset(df = ., it) %>%
                        lapply(function(x) {
                          x %>% mutate(ItemID = d)
                        })
                    })
    train <- lapply(split, function(x) {x$Training}) %>% bind_rows()
    test <- lapply(split, function(x) {x$Testing}) %>% bind_rows()
    
    # Deal with reso=0.1 
    d_subj <- model$item_spec %>% filter(Resolution == 0.1) %>% pull(ItemID)
    l <- lapply(list(train, test),
                function(x) {
                  x %>%
                    mutate(Resolution = case_when(ItemID %in% d_subj ~ 0.1,
                                                  TRUE ~ 1),
                           Score = round(Score / Resolution)) %>%
                    select(-Resolution)
                })

    # Calibration data
    if (model$include_calibration) {
      train_cal <- scorad %>% filter(Iteration <= it) %>%
        mutate(Resolution = case_when(Label %in% detail_POSCORAD("Subjective symptoms")$Label ~ 0.1,
                                      TRUE ~ 1),
               Score = round(Score / Resolution)) %>%
        select(-Resolution)
    } else {
      train_cal <- NULL
    }
    
    # Treatment data
    if (model$include_treatment) {
      train_treat <- treat %>% filter(Iteration <= it)
    } else {
      train_treat <- NULL
    }
    
    # Add recommendations input
    if (model$include_recommendations) {
      # Make recommendation at the last time of the training iteration (whether there is a training observation, or observed outcome)
      df_rec <- data.frame(Patient = pt,
                           Time = fc_it %>% filter(Iteration == it) %>% pull(LastTime)) %>%
        full_join(t_max, by = "Patient") %>%
        filter(Time <= LastTime) %>%
        select(-LastTime)
    } else {
      df_rec <- NULL
    }
    
    data_stan <- c(data_stan0,
                   prepare_standata(model,
                                    train = l[[1]], 
                                    test = l[[2]],
                                    cal = train_cal,
                                    treat = train_treat,
                                    rec = df_rec))
    
    id <- bind_rows(l[[1]], l[[2]], train_cal, train_treat, df_rec) %>%
      get_index()
    
    fit <- sampling(compiled_model,
                    data = data_stan,
                    pars = param,
                    control = list(adapt_delta = 0.9),
                    init = 0,
                    iter = n_it,
                    chains = n_chains,
                    refresh = 0)
    
    ## Performance of individual signs
    pred <- rstan::extract(fit, pars = "y_pred")[[1]]
    smp <- lapply(1:ncol(pred), function(i) {pred[, i]})
    perf <- test %>%
      mutate(lpd0 = extract_lpd(fit),
             Samples = smp)
    for (d in 1:model$D) {
      perf %>%
        filter(ItemID == d) %>%
        select(-ItemID) %>%
        mutate(Samples = map(Samples, ~(.x * model$item_spec$Resolution[d]))) %>%
        add_metrics2_d(support = seq(0, model$item_spec$Maximum[d], model$item_spec$Resolution[d])) %>%
        select(-lpd) %>%
        rename(lpd = lpd0) %>%
        saveRDS(file = here(file_dict[[model$item_spec$Name[d]]]$ValDir,
                                 paste0("val_", it, ".rds")))
    }
    
    ## Performance of aggregates
    pred_agg <- rstan::extract(fit, pars = "agg_rep")[[1]]
    agg_names <- gsub("weight_", "", colnames(data_stan$agg_weights))
    for (d in 1:length(agg_names)) {
      # Obtain test set for aggregate
      agg_dict <- detail_POSCORAD() %>%
        filter(Name == agg_names[d])
      test_agg <- POSCORAD %>%
        rename(Score = all_of(agg_dict$Label)) %>%
        select(Patient, Time, Score) %>%
        mutate(Iteration = get_fc_iteration(Time, t_horizon)) %>%
        split_fc_dataset(df = ., it)
      test_agg <- test_agg$Testing
      # Extract predictive samples
      id_test <- left_join(test_agg, id, by = c("Patient", "Time")) %>% pull(Index)
      smp_agg_d <- lapply(seq_along(id_test), function(i) {pred_agg[, id_test[i], d]})
      perf_agg <- test_agg %>%
        mutate(Samples = smp_agg_d) # replace by EczemaPred::samples_to_list(pred_agg[, id_test, d])
      if (agg_names[d] %in% c("SCORAD", "oSCORAD")) {
        perf_agg <- perf_agg %>%
          add_metrics2_c(., add_samples = 0:agg_dict$Maximum, bw = 0.5)
      } else {
        perf_agg <- perf_agg %>%
          add_metrics2_d(., support = seq(0, agg_dict$Maximum, agg_dict$Resolution))
      }
      # Save validation results (better to save in the loop in case something breaks)
      saveRDS(perf_agg, file = here(file_dict[[agg_names[d]]]$ValDir,
                                         paste0("val_", it, ".rds")))
    }
    
    ## Recommendations
    if (model$include_recommendations) {
      
      aggrec <- rstan::extract(fit, pars = "agg_rec")[[1]]
      yrec <- rstan::extract(fit, pars = "y_rec")[[1]]

      # Add severity item samples to pred_rec
      pred_rec <- df_rec
      for (d in 1:nrow(model$item_spec)) {
        tmp <- model$item_spec[d, ]
        pred_rec[[tmp$Label]] <- lapply(1:nrow(pred_rec),
                                        function(j) {
                                          yrec[, , j, tmp$ItemID]
                                        })
      }
      # Add aggregates samples to pred_rec
      for (d in seq_along(agg_names)) {
        pred_rec[[detail_POSCORAD(agg_names[d])$Label]] <- lapply(1:nrow(pred_rec),
                                                                  function(j) {
                                                                    aggrec[, , j, d]
                                                                  })
      }
      # Add p_treat to pred_rec
      df_rec <- left_join(df_rec, id, by = c("Patient", "Time"))
      ptreat <- rstan::extract(fit, pars = "p_treat")[[1]]
      ptreat <- ptreat[, df_rec[["Index"]], ]
      for (i in seq_along(model$treatment_names)) {
        pred_rec[[paste0(model$treatment_names[i], "_post")]] <- lapply(1:nrow(pred_rec),
                                                                  function(j) {
                                                                    ptreat[, j, i]
                                                                  })
      }
      
      # Save recommendation results
      saveRDS(list(Predictions = pred_rec, Actions = model$actions),
              file = here(rec_files$RecDir, paste0("rec_", it, ".rds")))
    }
    
    ####
    
    duration <- Sys.time() - duration
    cat(glue::glue("Ending iteration {it} after {round(duration, 1)} {units(duration)}"), sep = "\n")
    
    # Return
    NULL
  }
  stopCluster(cl)
  
  # Recombine validation results
  for (j in 1:length(file_dict)) {
    recombine_results(dir_name = file_dict[[j]]$ValDir,
                      output_file = file_dict[[j]]$Val,
                      expected_number_of_files = length(train_it))
  }
  
  # Recombine recommendation results
  if (model$include_recommendations) {
    
    # Check actions dataframes
    files <- list.files(rec_files$RecDir, full.names = TRUE)
    list_actions <- lapply(files,
                           function(x) {
                             tmp <- readRDS(x)
                             return(tmp[["Actions"]])
                           })
    all_same_actions <- all(vapply(list_actions, function(x) {all.equal(x, list_actions[[1]])}, logical(1)))
    if (!all_same_actions) {
      warning("The actions dataframes are not the same across iterations.")
    }
    
    recombine_results(dir_name = rec_files$RecDir,
                      output_file = rec_files$RecFile,
                      reading_function = function(x) {readRDS(x)[["Predictions"]]},
                      expected_number_of_files = length(train_it))
    res_rec <- readRDS(rec_files$RecFile)
    saveRDS(list(Predictions = res_rec, Actions = list_actions[[1]]),
            file = rec_files$RecFile)
    
  }
  
}
