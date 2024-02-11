# Notes -------------------------------------------------------------------

# Plot performance at a given iteration (as a pointrange) for all severity items and aggregates (x axis), models (colour)

# Initialisation ----------------------------------------------------------

rm(list = ls()) # Clear Workspace (better to restart the session)

source(here::here("analysis", "00_init.R"))

#### OPTIONS
metric <- "lpd" # for comparing learning curves, not paired comparisons
t_horizon <- 4 # horizon that was used for forward chaining
it <- 16 # Iteration to plot the performance for
max_horizon <- 14 # restrict prediction horizon
pred_horizon <- 4 # prediction horizon
####

metric <- match.arg(metric, c("lpd", "RPS"))
stopifnot(t_horizon > 0,
          max_horizon > 0)

dataset <- "PFDC"

item_dict <- detail_POSCORAD() %>%
  mutate(Component = case_when(Name == "extent" ~ "Extent",
                               Name %in% c("itching", "sleep") ~ "Subjective symptoms",
                               Name %in% detail_POSCORAD("Intensity signs")$Name ~ "Intensity signs",
                               TRUE ~ "Aggregates"))

list_models <- available_models() %>%
  select(Model) %>%
  mutate(ModelGroup = "EczemaPred")
list_models <- tibble(Model = c("uniform", "historical"),
                      ModelGroup = "Reference") %>%
  bind_rows(list_models)
mdl_names <- list_models[["Model"]]
list_models <- list_models %>%
  expand_grid(item_dict) %>%
  rename(Item = Name)

list_models[["File"]] <- vapply(1:nrow(list_models),
                                function(i) {
                                  get_results_files(outcome = list_models$Item[i],
                                                    model = list_models$Model[i],
                                                    dataset = dataset,
                                                    val_horizon = t_horizon,
                                                    root_dir = here())$Val
                                },
                                character(1))
stopifnot(all(file.exists(list_models$File)))

# Processing --------------------------------------------------------------

fc_it <- load_PFDC()$POSCORAD %>%
  rename(Time = Day) %>%
  detail_fc_training(df = ., horizon = t_horizon)

perf <- lapply(1:nrow(list_models),
               function(i) {
                 readRDS(list_models$File[i]) %>%
                   filter(Horizon <= max_horizon,
                          Iteration > 0 | list_models$Model[i] != "RW") %>%
                   estimate_performance(metric, ., fc_it, adjust_horizon = !(list_models$Model[i] %in% c("historical", "uniform"))) %>%
                   bind_cols(., list_models[i, ])
               }) %>%
  bind_rows()

# Performance at given iteration --------------------------------------------------------------------

pal <- rev(cbbPalette[1:length(unique(list_models[["Model"]]))])

tmp <- perf %>%
  filter(Variable == "Fit",
         Iteration == it,
         Horizon == pred_horizon) %>%
  mutate(ModelGroup = factor(ModelGroup, levels = rev(c("EczemaPred", "Reference"))),
         Model = factor(Model, levels = mdl_names))

brk <- c(.01, .1, .25, .5, 1)

## Performance items

item_names <- detail_POSCORAD("Items")$Name

p1 <- tmp %>%
  filter(Item %in% item_names) %>%
  ggplot(aes(x = Item, y = Mean, ymin = Mean - SE, ymax = Mean + SE, colour = Model, shape = ModelGroup)) +
  facet_grid(rows = vars(Component), space = "free", scale = "free") +
  geom_pointrange(position = position_dodge(width = .66), size = 1, fill = "white") +
  scale_colour_manual(values = pal) +
  scale_y_continuous(breaks = log(brk), labels = paste0("log(", brk, ")")) +
  scale_shape_manual(values = c(16, 21)) +
  coord_flip() +
  labs(x = "", y = metric, colour = "", shape = "") +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = .5, hjust = .5))
p1

## Performance for aggregate scores

p2 <- tmp %>%
  filter(Item %in% c("oSCORAD", "SCORAD")) %>%
  ggplot(aes(x = Item, y = Mean, ymin = Mean - SE, ymax = Mean + SE, colour = Model, shape = ModelGroup)) +
  facet_grid(rows = vars(Component), space = "free", scale = "free") +
  geom_pointrange(position = position_dodge(width = .66), size = 1, fill = "white") +
  scale_colour_manual(values = pal) +
  scale_shape_manual(values = c(16, 21)) +
  coord_flip() +
  labs(x = "", y = metric, colour = "", shape = "") +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = .5, hjust = .5))

## Custom legend

legend1 <-  tmp %>%
  filter(Item %in% item_names) %>%
  filter(ModelGroup == "EczemaPred") %>%
  ggplot(aes(x = Item, y = Mean, ymin = Mean - SE, ymax = Mean + SE, colour = Model)) +
  geom_pointrange(position = position_dodge(width = .66), size = 1, fill = "white", shape = 21) +
  scale_colour_manual(values = pal[-(1:2)], name = "EczemaPred models") +
  labs(x = "", y = metric, colour = "") +
  theme(legend.position = "top") +
  guides(colour = guide_legend(title.position = "top"))
legend1 <- get_legend(legend1)

legend2 <- tmp %>%
  filter(Item %in% item_names) %>%
  filter(ModelGroup == "Reference") %>%
  ggplot(aes(x = Item, y = Mean, ymin = Mean - SE, ymax = Mean + SE, colour = Model)) +
  geom_pointrange(position = position_dodge(width = .66), size = 1, shape = 16) +
  scale_colour_manual(values = pal[1:2], name = "Reference models") +
  labs(x = "", y = metric, colour = "") +
  theme(legend.position = "top") +
  guides(colour = guide_legend(title.position = "top"))
legend2 <- get_legend(legend2)

legend <- plot_grid(NULL, legend1, NULL, legend2, NULL, nrow = 1)

## Combine plots

plot_grid(legend,
          plot_grid(p1, p2,
                    labels = "AUTO",
                    ncol = 1, rel_heights = c(9, 2), align = "v"),
          ncol = 1,
          rel_heights = c(1, 11))

# ggsave(here("results", "performance.jpg"), width = 13, height = 8, units = "cm", dpi = 300, scale = 3.5, bg = "white")
