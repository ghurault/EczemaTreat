# Processing power prior -------------------------------------------------------------

#' Extract power prior from historical posterior summary statistics and put in nice format for univariate models
#'
#' @param par_historical Dataframe of summary statistics of population parameters (for only one model)
#'
#' @return Named list (prefixed by "historical_") of parameters in par_historical:
#' - the first element/column correspond to the mean
#' - the second element/column corresponds to the standard deviation.
extract_powerprior_uni <- function(par_historical) {
  # Can be useful to integrate to EczemaPred
  
  stopifnot(is.data.frame(par_historical),
            all(c("Variable", "Mean", "sd") %in% colnames(par_historical)))
  
  par_names <- unique(par_historical[["Variable"]])
  power_prior <- lapply(par_names,
                        function(x) {
                          tmp <- par_historical %>%
                            filter(Variable %in% x) %>%
                            arrange(Index)
                          out <- cbind(tmp[["Mean"]], tmp[["sd"]])
                          if (nrow(tmp) == 1) {
                            out <- c(out)
                          }
                          return(out)
                        })
  names(power_prior) <- paste0(par_names)
  
  stopifnot(all(do.call(rbind, power_prior)[, 2] >= 0))
  
  return(power_prior)
  
}

#' Extract and process power prior for multivariate models
#'
#' @param par_historical Dataframe containing posterior summary statistics for all items.
#' It should notably associate each item to a unique Item ID and Distribution ID.
#'
#' @return List containing power prior to add to Stan data input
process_powerprior_FullModel <- function(par_historical) {
  
  stopifnot(is.data.frame(par_historical),
            all(c("Distribution", "ItemID") %in% colnames(par_historical)))
  # Other columns are checked in extract_powerprior_uni
  
  comp1 <- par_historical %>% filter(Distribution == 1) %>% distinct(ItemID) %>% pull()
  D1 <-  length(comp1)
  comp2 <- par_historical %>% filter(Distribution == 2) %>% distinct(ItemID) %>% pull()
  D2 <-  length(comp2)
  D <- D1 + D2
  
  # Extract power prior for each item
  lpp <- lapply(1:D,
                function(d) {
                  par_historical %>%
                    filter(ItemID == d) %>%
                    extract_powerprior_uni()
                })
  
  # Everything but delta
  power_prior <- lapply(lpp,
                        function(x) {
                          x[["delta"]] <- NULL
                          return(x)
                        }) %>%
    merge_lists(along = 0)
  
  # Deal with delta
  tmp <- list(lpp[comp1], lpp[comp2])
  pp_delta <- lapply(1:length(tmp),
                     function(i) {
                       lapply(tmp[[i]],
                              function(x) {
                                out <- x[c("delta")]
                                names(out) <- paste0(names(out), i)
                                return(out)
                              }) %>%
                         merge_lists(along = 0)
                     }) %>%
    unlist(recursive = FALSE)
  
  power_prior <- c(pp_delta, power_prior)
  
  return(power_prior)
}

# Power prior influence ---------------------------------------------------

#' Associate forgetting factor `a0` to the relative importance `rho` of new dataset in the final posterior
#'
#' @param rho Relative importance of new dataset in the final posterior
#' @param n_new Number of observation in the new dataset.
#' Default to number of observations in PFDC dataset.
#' @param n_historical Number of observations in the historical dataset.
#' Default to the number of observations in the Derexyl dataset (=9408).
#' This value can be obtained with `nrow(load_dataset("Derexyl"))` 
#' in [EczemaPredPOSCORAD](https://github.com/ghurault/EczemaPredPOSCORAD).
#'
#' @return Vector
rho_to_a0 <- function(rho,
                      n_new = nrow(load_PFDC()$POSCORAD),
                      n_historical = 9408) {
  n_new / n_historical * (1 - rho) / rho
}

#' Associate the relative importance `rho` of new dataset in the final posterior to the forgetting factor `a0`
#'
#' @param a0 Forgetting factor
#' @param n_new Number of observation in the new dataset.
#' Default to number of observations in PFDC dataset.
#' @param n_historical Number of observations in the historical dataset.
#' Default to the number of observations in the Derexyl dataset (=9408).
#' This value can be obtained with `nrow(load_dataset("Derexyl"))` 
#' in [EczemaPredPOSCORAD](https://github.com/ghurault/EczemaPredPOSCORAD).
#'
#' @return Vector
a0_to_rho <- function(a0,
                      n_new = nrow(load_PFDC()$POSCORAD),
                      n_historical = 9408) {
  n_new / (n_new + a0 * n_historical)
}

#' Plot influence of forgetting factor when using a power prior
#' 
#' Plot the relative importance `rho` of new dataset in the final posterior as a function of the forgetting factor `a0`.
#'
#' @param n_new Number of observation in the new dataset.
#' Default to number of observations in PFDC dataset.
#' @param n_historical Number of observations in the historical dataset.
#' Default to the number of observations in the Derexyl dataset (=9408).
#' This value can be obtained with `nrow(load_dataset("Derexyl"))` 
#' in [EczemaPredPOSCORAD](https://github.com/ghurault/EczemaPredPOSCORAD).
#'
#' @return Ggplot
plot_powerprior_influence <- function(n_new = nrow(load_PFDC()$POSCORAD),
                                      n_historical = 9408) {
  tibble(a0 = seq(0, 1, .01)) %>%
    mutate(rho = a0_to_rho(a0)) %>%
    ggplot(aes(x = a0, y = rho)) +
    geom_line() +
    coord_cartesian(xlim = c(0,1 ), ylim = c(0, 1), expand = FALSE) +
    theme_bw(base_size = 15)
}
