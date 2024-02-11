# Notes -------------------------------------------------------------------

# Save posterior summary statistics of the main parameters, for all items, for the models trained with the Derexyl dataset.

# This script assume posterior summary statistics of the OrderedRW (v2) model fitted with the Derexyl dataset
# was saved at the location given by `get_results_files`.
# The script to fit the OrderedRW model can be found in [EczemaPredPOSCORAD](https://github.com/ghurault/EczemaPredPOSCORAD).

# Initialisation ----------------------------------------------------------

rm(list = ls()) # Clear Workspace (better to restart the session)

source(here::here(here::here("analysis", "00_init.R")))

intensity_signs <- detail_POSCORAD("Intensity signs")$Name
subjective_symptoms <- detail_POSCORAD("Subjective symptoms")$Name

####
datasets <- c("Derexyl")
####

param <- list_parameters("OrderedRW")

list_models <- data.frame(Item = detail_POSCORAD("Items")$Name,
                          Model = "OrderedRW") %>%
  mutate(PopulationParameter = list(param$Population),
         PatientParameter = list(param$Patient))
l_prior <- list_models %>%
  mutate(Distribution = "Prior",
         File = map2(Item, Model, ~get_results_files(outcome = .x, model = .y, root_dir = here())$PriorPar) %>% unlist())

l_post <- list_models %>%
  expand_grid(., Dataset = datasets) %>%
  mutate(Distribution = paste0("Posterior - ", Dataset),
         File = pmap(list(Item, Model, Dataset),
                     function(x, y, z) {
                       get_results_files(outcome = x, model = y, dataset = z, root_dir = here())$FitPar
                     }) %>%
           unlist()) %>%
  select(-Dataset)

list_models <- bind_rows(l_prior, l_post)

stopifnot(all(file.exists(list_models$File)))

# Combine estimates -------------------------------------------------------

par_POSCORAD <- lapply(unique(list_models[["Item"]]),
              function(item) {
                tmp <- list_models %>% filter(Item == item)

                out <- lapply(1:nrow(tmp),
                              function(i) {
                                readRDS(tmp$File[i]) %>%
                                  mutate(Distribution = tmp$Distribution[i],
                                         Model = tmp$Model[i]) %>%
                                  filter(Variable %in% tmp$PopulationParameter[i][[1]] |
                                           (Distribution == "Prior" & Variable %in% tmp$PatientParameter[i][[1]] & Index == 1))
                              }) %>%
                  bind_rows() %>%
                  select(-Patient, -Time) %>%
                  mutate(Distribution = forcats::fct_relevel(factor(Distribution), "Prior"),
                         Item = item)

                return(out)
              }) %>%
  bind_rows()

saveRDS(par_POSCORAD, file = here("data", "par_POSCORAD.rds"))

# Compare fit -------------------------------------------------------------

# par_POSCORAD <- EczemaPredPOSCORAD::par_POSCORAD

bind_rows(par_POSCORAD %>%
            filter(!is.na(Index)) %>%
            mutate(Variable = paste0(Variable, "[", Index, "]")),
          par_POSCORAD %>%
            filter(is.na(Index))) %>%
  select(-Index) %>%
  ggplot(aes(x = Variable, y = Mean, ymin = `5%`, ymax = `95%`, colour = Distribution)) +
  facet_wrap(vars(Item), scales = "free") +
  geom_pointrange(position = position_dodge(width = .25)) +
  coord_flip() +
  scale_colour_manual(values = cbbPalette) +
  labs(x = "", y = "Estimate", colour = "") +
  theme_bw(base_size = 15) +
  theme(legend.position = "top")
# ggsave(here("results", "compare_fit.jpg"), width = 13, height = 8, units = "cm", dpi = 300, scale = 3)
