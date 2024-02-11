# Notes -------------------------------------------------------------------

# Plot calibration estimates and plot SCORAD

# Initialisation ----------------------------------------------------------

rm(list = ls()) # Clear Workspace (better to restart the session)

source(here::here("analysis", "00_init.R"))

#### OPTIONS
pid <- 3
####

model <- ScoradPred(independent_items = FALSE,
                    a0 = .04,
                    include_trend = FALSE,
                    include_calibration = TRUE,
                    include_treatment = TRUE,
                    treatment_names = c("localTreatment", "emollientCream"),
                    include_recommendations = FALSE)

file_dict <- get_results_files(outcome = "SCORAD",
                               model = model$name,
                               dataset = "PFDC",
                               root_dir = here())

# Load data ---------------------------------------------------------------

l <- load_PFDC()

POSCORAD <- l$POSCORAD %>%
  rename(Time = Day)
df <- POSCORAD %>%
  select(one_of("Patient", "Time", model$item_spec$Label)) %>%
  pivot_longer(cols = all_of(model$item_spec$Label), names_to = "Label", values_to = "Score") %>%
  drop_na()  %>%
  left_join(model$item_spec[, c("Label", "ItemID")], by = c("Label"))

# Prepare SCORAD calibration data
if (model$include_calibration) {
  cal <- scorad <- l$SCORAD %>%
    rename(Time = Day) %>%
    select(one_of("Patient", "Time", model$item_spec$Label)) %>%
    pivot_longer(cols = all_of(model$item_spec$Label), names_to = "Label", values_to = "Score") %>%
    drop_na()  %>%
    left_join(model$item_spec[, c("Label", "ItemID")], by = c("Label")) 
} else {
  cal <- NULL
}

# Prepare treatment data
treatment_lbl <- paste0(model$treatment_names, "WithinThePast2Days")
if (model$include_treatment) {
  treat <- POSCORAD %>%
    select(all_of(c("Patient", "Time", treatment_lbl))) %>%
    pivot_longer(cols = all_of(treatment_lbl), names_to = "Treatment", values_to = "UsageWithinThePast2Days") %>%
    mutate(Treatment = vapply(Treatment, function(x) {which(x == treatment_lbl)}, numeric(1)) %>% as.numeric()) %>%
    drop_na()
} else {
  treat <- NULL
}

# NB: assume no recommendation (at least outside time-series)

pt <- unique(df[["Patient"]])

id <- get_index(bind_rows(df, cal, treat))
df <- left_join(df, id, by = c("Patient", "Time"))

# Results
fit <- readRDS(file_dict$Fit)
par <- readRDS(file_dict$FitPar)

# Correlation plot --------------------------------------------------------

x <- "Omega"

omg <- rstan::extract(fit, pars = x)[[1]]

tmp <- list(Mean = apply(omg, c(2, 3), mean),
            SD = apply(omg, c(2, 3), sd),
            Lower = apply(omg, c(2, 3), function(x) {quantile(x, probs = .05)}),
            Upper = apply(omg, c(2, 3), function(x) {quantile(x, probs = .95)}),
            pval = apply(omg, c(2, 3), function(x) {empirical_pval(x, 0)}))
tmp <- lapply(tmp,
              function(x) {
                colnames(x) <- model$item_spec$Name
                rownames(x) <- model$item_spec$Name
                return(x)
              })

jpeg(here("results", paste0(x, "_", model$name, ".jpeg")),
     width = 20, height = 20, units = "cm", res = 300, quality = 95, pointsize = 11)
corrplot::corrplot.mixed(tmp$Mean, lower = "number", upper = "ellipse")
dev.off()

# Combine with power prior plot in `plot_powerprior.R`

# Calibration plot --------------------------------------------------------------------

# Estimates
p1_cal <- par %>%
  filter(Variable == "bias0") %>%
  rename(ItemID = Index) %>%
  left_join(model$item_spec, by = "ItemID") %>%
  filter(!(Name %in% c("sleep", "itching"))) %>%
  mutate(Name = factor(Name),
         Name = factor(Name, levels = rev(levels(Name)))) %>%
  ggplot(aes(x = Name, y = Mean, ymin = `5%`, ymax = `95%`)) +
  facet_grid(rows = vars(Component), scales = "free", space = "free") +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  scale_y_continuous(limits = c(-.5, .5),
                     breaks = c(-.5, -.25, 0, .25, 0.5),
                     labels = c("-0.5\nPatient scores\nhigher than clinician",
                                -0.25, 0, 0.25,
                                "0.5\nClinician scores\n higher than patient")) +
  labs(x = "",
       y = "Initial bias (normalised)")

# Otherwise, post-process figures to give the interpretation of the direction of the effect
# ("patient scores higher than clinician" vs "clinician scores higher than patient")

aggcal <- rstan::extract(fit, pars = "agg_cal_rep")[[1]]
aggcal <- aggcal[, , 4] # SCORAD

### Plot observed PO-SCORAD and inferred SCORAD as a fanchart
tmp <- POSCORAD %>%
  filter(Patient == pid)
p2_cal <- plot_post_traj_fanchart(aggcal,
                                  id = id,
                                  patient_id = pid,
                                  legend_fill = "discrete",
                                  CI_level = seq(0.1, 0.9, 0.2),
                                  max_score = 60) +
  add_broken_pointline(tmp, aes_x = "Time", aes_y = "SCORAD", colour = "Observed\nPO-SCORAD") +
  scale_colour_manual(values = c("Observed\nPO-SCORAD" = "black")) +
  labs(fill = "Inferred\nSCORAD\nprobabilities", colour = "") +
  theme(legend.position = c(.9, .8),
        legend.title = element_text(size = 11),
        legend.spacing.y = unit(0, 'cm'))

plot_grid(p1_cal, p2_cal, nrow = 1, labels = "AUTO")

if (FALSE) {
  ggsave(here("results", "plot_calibration.jpg"),
         width = 18, height = 7, units = "cm", dpi = 300, scale = 2.5)
}

# Treatment ---------------------------------------------------------------

p_treat <- extract_par_indexes(par, var_name = "ATE", dim_names = c("ItemID", "Treatment")) %>%
  filter(Variable == "ATE") %>%
  mutate(Treatment = model$treatment_names[Treatment]) %>%
  left_join(model$item_spec, by = "ItemID") %>%
  mutate(Treatment = recode(Treatment,
                            emollientCream = "Emollient Cream",
                            localTreatment = "Topical Corticosteroids"),
         Component = gsub(" ", "\n", Component),
         Name = factor(Name),
         Name = factor(Name, levels = rev(levels(Name)))) %>%
  ggplot(aes(x = Name, y = Mean, ymin = `5%`, ymax = `95%`, colour = Treatment)) +
  facet_grid(rows = vars(Component), scale = "free", space = "free") +
  geom_pointrange(position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  scale_y_continuous(limits = c(-.05, .05),
                     breaks = c(-.05, -.025, 0, .025, 0.05),
                     labels = c("-0.05\nTreatment\nreduces severity", -0.025, 0, 0.025, "0.5\nTreatment\nincreases severity")) +
  scale_colour_manual(values = cbbPalette[c(2, 1)]) +
  labs(x = "", y = "Treatment effect (normalised)", colour = "") +
  theme(legend.position = "top")
p_treat
# ggsave(here("results", "treatment_effects.jpg"), width = 13, height = 8, units = "cm", scale = 2)

# Combine with recommendation plot
p_rec <- readRDS(here("results", "subplot_recommendation.rds")) +
  labs(title = "")
plot_grid(p_treat, p_rec, nrow = 1, labels = "AUTO")

if (FALSE) {
  ggsave(here("results", "plot_treatment.jpg"),
         width = 10, height = 5, units = "cm", dpi = 300, scale = 3.5)
}
