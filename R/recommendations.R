# Recommendations ---------------------------------------------------------

#' Action dataframe
#'
#' @param names Vector of action names
#'
#' @return Dataframe with column from `names` and a column `Action`
get_actions <- function(names) {
  stopifnot(is.vector(names, mode = "character"))
  out <- lapply(seq_along(names), function(x) {0:1}) %>%
    expand.grid()
  colnames(out) <- names
  out[["Action"]] <- 1:nrow(out)
  return(out)
}

#' Dictionary of recommendation result files
#'
#' @param outcome Score to predict
#' @param model Model name
#' @param dataset Dataset name
#' @param val_horizon Validation horizon
#' @param root_dir Root directory. Default to working directory.
#'
#' @return Named list with elements "RecDir" and "RecFile"
get_recommendation_files <- function(outcome, model, dataset, val_horizon, root_dir = ".") {
  check_input <- get_results_files(outcome = outcome, model = model, dataset = dataset, val_horizon = val_horizon)
  
  out <- list(RecDir = file.path(root_dir, "results", glue::glue("rec{val_horizon}_{outcome}_{model}_{dataset}")))
  out$RecFile <- paste0(out$RecDir, ".rds")
  
  return(out)
}
