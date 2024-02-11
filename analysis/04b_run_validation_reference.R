# Notes -------------------------------------------------------------------

# Run validation for univariate reference models (uniform, historical forecast)

# Initialisation ----------------------------------------------------------

rm(list = ls()) # Clear Workspace (better to restart the session)

set.seed(2021) # Reproducibility (Stan use a different seed)

source(here::here("analysis", "00_init.R"))
library(foreach)
library(doParallel)

dataset <- "PFDC"

#### OPTIONS
score <- "SCORAD"
mdl_name <- "historical"

run <- FALSE
t_horizon <- 4
n_chains <- 4
n_it <- 2000
n_cluster <- 2 # floor((parallel::detectCores() - 2) / n_chains)
####

item_dict <- detail_POSCORAD()

score <- match.arg(score, item_dict[["Name"]])
mdl_name <- match.arg(mdl_name, c("uniform", "historical"))

item_dict <- item_dict %>% filter(Name == score)
item_lbl <- as.character(item_dict[["Label"]])
max_score <- item_dict[["Maximum"]]
reso <- item_dict[["Resolution"]]
M <- round(max_score / reso)

is_continuous <- (score %in% c("SCORAD", "oSCORAD"))

## Files
file_dict <- get_results_files(outcome = score,
                               model = mdl_name,
                               dataset = dataset,
                               val_horizon = t_horizon,
                               root_dir = here())

# Data ---------------------------------------------------------------------

POSCORAD <- load_PFDC()$POSCORAD

# Subset dataset
df <- POSCORAD %>%
  rename(Time = Day, Score = all_of(item_lbl)) %>%
  select(Patient, Time, Score) %>%
  drop_na()

pt <- unique(df[["Patient"]])

# Forward chaining --------------------------------------------------------

df <- df %>% mutate(Iteration = get_fc_iteration(Time, t_horizon))
train_it <- get_fc_training_iteration(df[["Iteration"]])

if (run) {
  
  cl <- makeCluster(n_cluster, outfile =  "")
  registerDoParallel(cl)
  
  dir.create(file_dict$ValDir)
  
  out <- foreach(i = rev(seq_along(train_it))) %dopar% {
    it <- train_it[i]
    
    # Need to reload functions and libraries
    source(here::here("analysis", "00_init.R"))
    
    duration <- Sys.time()
    cat(glue::glue("Starting iteration {it}"), sep = "\n")
    
    ####
    
    split <- split_fc_dataset(df, it)
    train <- split$Training
    test <- split$Testing
    
    # Uniform forecast
    if (mdl_name == "uniform" && !is_continuous) {
      perf <- test %>%
        mutate(Score = round(Score / reso)) %>%
        add_uniform_pred(test = .,
                         max_score = M,
                         discrete = TRUE,
                         include_samples = FALSE) %>%
        mutate(Score = Score * reso)
    }
    if (mdl_name == "uniform" && is_continuous) {
      perf <- test %>%
        add_uniform_pred(test = .,
                         max_score = max_score,
                         discrete = FALSE,
                         include_samples = TRUE,
                         n_samples = 2 * max_score)
    }
    
    # Historical forecast
    if (mdl_name == "historical" && !is_continuous) {
      perf <- test %>%
        mutate(Score = round(Score / reso)) %>%
        add_historical_pred(test = .,
                            train = mutate(train, Score = round(Score / reso)),
                            max_score = M,
                            discrete = TRUE,
                            add_uniform = TRUE,
                            include_samples = FALSE) %>%
        mutate(Score = Score * reso)
    }
    if (mdl_name == "historical" && is_continuous) {
      perf <- test %>%
        add_historical_pred(test = .,
                            train = train,
                            max_score = max_score,
                            discrete = FALSE,
                            add_uniform = TRUE,
                            include_samples = TRUE)
    }

    perf <- perf %>%
      select(-LastTime, -LastScore)
    
    # Save results (better to save in the loop in case something breaks)
    saveRDS(perf, file = here(file_dict$ValDir, paste0("val_", it, ".rds")))
    
    ####
    
    duration <- Sys.time() - duration
    cat(glue::glue("Ending iteration {it} after {round(duration, 1)} {units(duration)}"), sep = "\n")
    
    # Return
    NULL
  }
  stopCluster(cl)

  # Recombine results
  files <- list.files(file_dict$ValDir, full.names = TRUE)
  if (length(files) < length(train_it)) {
    warning(glue::glue("Number of files (={length(files)}) less than the number of unique iterations (={length(train_it)}). 
                       Some runs may have failed."))
  }
  res <- lapply(files,
                function(f) {
                  readRDS(f)
                }) %>%
    bind_rows()
  saveRDS(res, file = file_dict$Val)
  
}
