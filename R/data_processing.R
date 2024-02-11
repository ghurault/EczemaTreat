# Fixing treatment in PFDC dataset ----------------------------------------

#' Impute patterns "0-NA-0", "0-NA-NA-0" by 0 in "treatment usage within the past two days"
#'
#' @param x Vector of binary values corresponding to the time-series of "treatment usage within the past two days"
#'
#' @return x after imputation
impute_tw2_missing <- function(x) {

  tibble(x = x) %>%
    mutate(Missing = is.na(x),
           is0 = (x == 0),
           Impute0 = (Missing & lag(is0) & lead(is0)) |
             (Missing & lag(is0) & lead(Missing) & lead(is0, 2)) |
             (Missing & lag(Missing) & lag(is0, 2) & lead(is0)),
           Impute0 = replace_na(Impute0, FALSE),
           x = replace(x, Impute0, 0)) %>%
    pull(x)
  
}

#' Find pattern in string and return logical vector of matches
#'
#' @param x String
#' @param pattern Pattern to match in `x`
#'
#' @return Return logical vector of length `nchar(x)` indicating the matches
#' 
#' @examples
#' find_pattern(x = "0010110100", pattern = "010")
find_pattern <- function(x, pattern) {

  stopifnot(is_scalar_character(x),
            is_scalar_character(pattern))
  
  matches <- gregexpr(pattern, x)[[1]]
  
  out <- rep(FALSE, nchar(x))
  for (i in seq_along(matches)) {
    if (matches[i] > 0) {
      out[matches[i]:(matches[i] + attr(matches, "match.length")[i] - 1)] <- TRUE
    }
  }
  
  return(out)
  
}

#' Replace contradictions in "treatment usage within the past two days" time-series by missing values
#'
#' @param x Vector of binary values corresponding to the time-series of treatment usage within the past two days
#'
#' @return `x` after replacing contradictions with NA
replace_tw2_contradictions <- function(x) {

  txt <- replace(x, is.na(x), "n") %>%
    paste0(collapse = "")
  contr <- find_pattern(txt, "010") | find_pattern(txt, "0((1[1n])|[1n]1)0")
  x <- replace(x, contr, NA)
  
  return(x)
  
}

#' Fix treatment time-series in PFDC dataset
#' 
#' For both localTreatment/emollientCream within the past two days:
#' - Impute missing values in patterns "0-NA-0" and "0-NA-NA-0" to 0
#' - Set illogical patterns "0-1-0" and "0-1-1-0" (with potentially one of the 1s missing) to missing
#' - Reimpute patterns "0-NA-0" and "0-NA-NA-0" to 0
#' 
#' On average, 10% of values are changed.
#' Note the output dataset as more rows than the input dataset (inferred treatment when originally absent).
#'
#' @param df POSCORAD PFDC dataset
#'
#' @return POSCORAD PFDC dataset with updated treatment columns
fix_treatment <- function(df) {
  
  lapply(unique(df[["Patient"]]),
         function(pid) {
           
           sub_df <- df %>%
             filter(Patient == pid)
           
           for (tname in c("localTreatmentWithinThePast2Days", "emollientCreamWithinThePast2Days")) {
             
             obs <- sub_df %>%
               rename(tw2 = all_of(tname)) %>%
               select(Patient, Day, tw2)
             
             mis <- tibble(Patient = pid,
                           Day = setdiff(1:max(obs[["Day"]]), obs[["Day"]]),
                           tw2 = NA)
             
             tmp <- bind_rows(obs, mis) %>%
               arrange(Day) %>%
               mutate(tw2 = tw2 %>%
                        impute_tw2_missing() %>%
                        replace_tw2_contradictions() %>%
                        impute_tw2_missing()) %>%
               drop_na()
             
             sub_df <- full_join(sub_df, tmp, by = c("Patient", "Day")) %>%
               select(!all_of(tname)) %>%
               rename_with(~replace(.x, .x == "tw2", tname))
             
           }
           
           return(sub_df)
           
         }) %>%
    bind_rows() %>%
    # Calculate date associated with new observations
    group_by(Patient) %>%
    mutate(Date = case_when(is.na(Date) ~ min(Date, na.rm = TRUE) + Day - 1,
                            TRUE ~ Date)) %>%
    ungroup()
  
}

# Processing PFDC dataset -------------------------------------------------

load_PFDC <- function() {
  # Load PFDC dataset
  # - Process treatment data
  # - Only consider patients with POSCORAD data
  # - Recompute day from by considering start date in POSCORAD and SCORAD datasets
  # - Remove patients with less than 5 observations
  # - Reindex patients between 1 and length(Patient)
  #
  # Args: none
  #
  # Returns:
  # Named list containing POSCORAD and SCORAD datasets
  
  poscorad <- TanakaData::POSCORAD_PFDC %>%
    fix_treatment() %>%
    arrange(Patient, Day)
  scorad <- TanakaData::SCORAD_PFDC
  
  poscorad <- factor_to_numeric(poscorad, "Patient")
  
  # Only consider patients with POSCORAD
  pt_poscorad <- unique(poscorad[["Patient"]])
  scorad <- scorad %>%
    filter(Patient %in% pt_poscorad)
  
  # Find start date
  start_date_scorad <- scorad %>%
    group_by(Patient) %>%
    summarise(StartDate = min(Date, na.rm = TRUE)) %>%
    ungroup()
  start_date_poscorad <- poscorad %>%
    group_by(Patient) %>%
    summarise(StartDate = min(Date, na.rm = TRUE)) %>%
    ungroup()
  start_date <- bind_rows(start_date_scorad, start_date_poscorad) %>%
    group_by(Patient) %>%
    summarise(StartDate = min(StartDate)) %>%
    ungroup()
  
  # Recompute Day in POSCORAD dataset
  poscorad <- full_join(poscorad, start_date, by = "Patient") %>%
    group_by(Patient) %>%
    mutate(Day = as.numeric(Date - StartDate + 1)) %>%
    ungroup()
  
  # Compute Day in SCORAD dataset
  scorad <- left_join(scorad, start_date, by = "Patient") %>%
    group_by(Patient) %>%
    mutate(Day = as.numeric(Date - StartDate + 1)) %>%
    ungroup()
  
  # Filter out patients and reindex them
  poscorad <- poscorad %>%
    group_by(Patient) %>%
    filter(n() >= 5) %>%
    mutate(NewID = cur_group_id()) %>%
    ungroup()
  patient_id <- poscorad %>%
    select(Patient, NewID) %>%
    distinct()
  scorad <- full_join(scorad, patient_id, by = "Patient") %>%
    mutate(Patient = NewID) %>%
    select(-NewID)
  poscorad <- poscorad %>%
    mutate(Patient = NewID) %>%
    select(-NewID)
  
  return(list(POSCORAD = poscorad, SCORAD = scorad))
}
