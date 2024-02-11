# Libraries and variables -------------------------------------------------

library(HuraultMisc) # Functions shared across projects
library(TanakaData)
library(EczemaPred)
library(EczemaPredPOSCORAD)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_bw(base_size = 15))
library(purrr)
library(magrittr)
library(here)
library(cowplot)

library(rstan)
rstan_options(auto_write = TRUE) # Save compiled model
options(mc.cores = parallel::detectCores()) # Parallel computing

lapply(list.files(here("R"), pattern = ".R$", full.names = TRUE),
       source)
