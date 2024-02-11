# Notes -------------------------------------------------------------------

# - Compare prior, historical posterior, posterior with/without power prior for a given parameter
# - For a given a0, how much does power prior influences posterior as a function of the number of data included
# The second plot is probably easier to understand.

# a0 indicates how informative the historical prior is,
# or how much information is shared between the historical posterior and the final posterior
# (very informative: weighted equally in the posterior; not informative: does not contribute to posterior)

# Initialisation ----------------------------------------------------------

rm(list = ls()) # Clear Workspace (better to restart the session)

source(here::here(here::here("analysis", "00_init.R")))

library(ggrepel)

#### OPTIONS
par_name <- "rho2" # check that processing of par2 works for something else than "rho2" (Index much match ItemID)
####

item_dict <- ScoradPred_items()

# Prior and posterior -----------------------------------------------------

# Prior and Posterior Derexyl
par1 <- readRDS(here("data", "par_POSCORAD.rds")) %>%
  filter(Variable == par_name)

# Posterior ScoradPred
par2 <- lapply(c("ScoradPred", "ScoradPred+h004"),
               function(m) {
                 get_results_files(outcome = "SCORAD", model = m, dataset = "PFDC", root_dir = here())$FitPar %>%
                   readRDS() %>%
                   mutate(Distribution = m)
               }) %>%
  bind_rows() %>%
  filter(Variable == par_name) %>%
  left_join(item_dict %>%
              select(Name, ItemID),
            by = c("Index" = "ItemID")) %>%
  rename(Item = Name)

par <- bind_rows(par1, par2) %>%
  mutate(Distribution = factor(Distribution,
                               levels = c("Prior", "Posterior - Derexyl", "ScoradPred", "ScoradPred+h004"),
                               labels = c("Prior", "Historical", "Without power prior", "With power prior")))

# Plot comparison estimate --------------------------------------------------------------------

par %>%
  select(-Index) %>%
  filter(Variable == par_name) %>%
  ggplot(aes(x = Item, y = Mean, ymin = `5%`, ymax = `95%`, colour = Distribution)) +
  geom_pointrange(position = position_dodge(width = .5)) +
  coord_flip() +
  scale_colour_manual(values = HuraultMisc::cbbPalette) +
  labs(x = "", y = par_name, colour = "") +
  theme(legend.position = "top")

# Plot relative importance ------------------------------------------------

poscorad <- load_PFDC()$POSCORAD %>%
  rename(Time = Day)

fc_it <- detail_fc_training(poscorad, horizon = 4)
id_xbrk2 <- vapply(seq(0, 1, length.out = 10), function(x) {which.min((x - fc_it$Proportion)^2)}, numeric(1))

a0_val <- c(.04,
            signif(rho_to_a0(.5), 2),
            .01, #signif(rho_to_a0(.9), 2),
            .5,
            1)

tmp <- expand_grid(a0 = a0_val,
                   N = 0:nrow(poscorad)) %>%
  mutate(rho = a0_to_rho(a0 = a0, n_new = N))

p <- tmp %>%
  group_by(a0) %>%
  mutate(Size = (a0 == .04)) %>%
  mutate(Label = ifelse(N == max(N), paste0("a0 = ", a0), NA)) %>%
  ungroup() %>%
  ggplot(aes(x = N, y = rho, colour = factor(a0, levels = a0_val), label = Label)) +
  geom_line(aes(size = Size)) +
  geom_label_repel(na.rm = TRUE) +
  labs(x = "Number of observations in training set",
       y = "Contribution to posterior (%)") +
  scale_x_continuous(expand = c(0, 0),
                     sec.axis = dup_axis(breaks =  fc_it$N[id_xbrk2],
                                         labels = fc_it$LastTime[id_xbrk2],
                                         name = "Training day")) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_colour_manual(values = HuraultMisc::cbbPalette[-c(5, 6)]) +
  scale_size_discrete(range = c(.5, 2)) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(legend.position = "none")
p
# ggsave(here("results", "powerprior.jpg"), width = 13, height = 8, units = "cm", dpi = 300, scale = 2)

if (FALSE) {
  # Combine power prior plot and correlogram for paper
  
  p_corr <- ggdraw() +
    draw_image(here("results", "Omega_ScoradPred+h004+corr+cal+treat.jpeg"),
               scale = 1)
  
  plot_grid(p, p_corr, nrow = 1, labels= "AUTO")
  ggsave(here("results", "powerprior_correlation.jpg"),
         width = 10, height = 5, units = "cm", dpi = 300, bg = "white", scale = 3.2)
  # Note that scale is for the first plot (not the image)
  
}
