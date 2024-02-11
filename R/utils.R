# Utility functions -----------------------------------------------------------

#' Full names of a two-dimensional parameter
#'
#' @param par_name Name of the parameter
#' @param n1 First dimension
#' @param n2 Second dimension
#'
#' @return Character vector
expand_2d_parname <- function(par_name, n1, n2) {
  expand_grid(i = 1:n1, j = 1:n2) %>%
    mutate(x = paste0(par_name, "[", i, ",", j, "]")) %>%
    pull(x)
}

#' Merge list of lists by names
#' 
#' Only elements with the same names are considered.
#' List's elements with no names in `ll` are discarded.
#' 
#' When `along = 0`, arrays are bind on a new dimension before the first:
#' for example, if the list all contains an element "a" which is a vector,
#' element "a" of the new list will be a matrix where rows will indicate the original list ID.
#'
#' @param ll List of lists
#' @param ... Arguments to pass to [abind::abind()]
#'
#' @return List
merge_lists <- function(ll, ...) {
  
  stopifnot(is.list(ll))
  
  lnames <- lapply(ll, names) %>% do.call(c, .) %>% unique()
  
  out <- lapply(lnames,
                function(nm) {
                  lapply(ll,
                         function(x) {
                           x[[nm]]
                         }) %>%
                    abind::abind(...)
                })
  names(out) <- lnames
  
  return(out)
}

#' Recombine results in different files
#' 
#' Read files corresponding to different dataframe, concatenate them and save it
#'
#' @param dir_name Directory containing intermediate files
#' @param output_file File to write the recombined results
#' @param reading_function Function used to read the files
#' @param expected_number_of_files (optional) expected number of files in `dir_name`
#'
#' @return NULL
recombine_results <- function(dir_name, output_file, reading_function = readRDS, expected_number_of_files = NULL) {
  
  stopifnot(is_scalar_character(dir_name),
            dir.exists(dir_name),
            is_scalar(output_file),
            is.function(reading_function))
  
  files <- list.files(dir_name, full.names = TRUE)
  
  if (!is.null(expected_number_of_files)) {
    stopifnot(is_scalar_wholenumber(expected_number_of_files))
    if (length(files) != expected_number_of_files) {
      warning(glue::glue("Number of files (={length(files)}) is different from the number of expected files (={expected_number_of_files})."))
    }
  }
  
  res <- lapply(files, reading_function) %>%
    bind_rows()
  saveRDS(res, file = output_file)
  
  return(NULL)
}

#' Extract index of `var_name` in `par`
#'
#' @param par Dataframe
#' @param var_name Character of the variable name to extract
#' @param dim_names Character vector of dimension names
#'
#' @return `par` with additional columns `dim_names`
extract_par_indexes <- function(par, var_name, dim_names) {
  
  stopifnot(is.data.frame(par),
            all(c("Variable") %in% colnames(par)),
            is_scalar_character(var_name),
            is.character(dim_names))
  
  rows_id <- grepl(paste0("^", var_name, "\\["), par[["Variable"]])
  sub_par <- filter(par, rows_id)
  
  id_sub <- HuraultMisc::extract_index_nd(sub_par[["Variable"]], dim_names = dim_names)
  
  sub_par[, c("Variable", "Index", dim_names)] <- NULL
  sub_par <- bind_cols(sub_par, id_sub)
  par <- par %>%
    filter(!rows_id) %>%
    bind_rows(sub_par)
  
  return(par)
}

# Available models --------------------------------------------------------

#' List models that are investigated
#'
#' @return Dataframe
available_models <- function() {
  data.frame(
    a0 = c(0, 0.04, 0.04, 0.04, 0.04, 0.04),
    independent_items = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
    include_calibration = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
    include_treatment = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
    include_trend = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)
  ) %>%
    mutate(Score = "SCORAD",
           Dataset = "PFDC",
           Args = pmap(list(a0 = a0,
                            independent_items = independent_items,
                            include_calibration = include_calibration,
                            include_treatment = include_treatment,
                            include_trend = include_trend),
                       list),
           Args = map(Args, ~c(.x, list(treatment_names = c("localTreatment", "emollientCream")))),
           Model = map(Args, ~do.call(ScoradPred, .x)$name) %>% unlist())
}
