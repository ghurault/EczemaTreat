# Notes -------------------------------------------------------------------

# Plot PO-SCORAD, SCORAD (breakdown) and treatment

# Initialisation ----------------------------------------------------------

rm(list = ls()) # Clear Workspace (better to restart the session)

source(here::here("analysis", "00_init.R"))

#### OPTIONS
pid <- 3
####

l <- load_PFDC()

pt <- unique(l$POSCORAD[["Patient"]])
stopifnot(pid %in% pt)

poscorad <- l$POSCORAD %>% filter(Patient == pid)
scorad <- l$SCORAD %>% filter(Patient == pid)

palette <- c("PO-SCORAD" = "#000000", "SCORAD" = "#E69F00")

# Plot breakdown --------------------------------------------------------------

item_dict <- detail_POSCORAD("Items")
subj_symp <- detail_POSCORAD("Subjective symptoms")$Label
intensity_signs <- detail_POSCORAD("Intensity signs")$Label

pl <- lapply(1:nrow(item_dict),
             function(i) {
               score_lbl <- item_dict$Label[i]
               score_lbl2 <- case_when(score_lbl == "Itching VAS" ~ "Itching",
                                       score_lbl == "Sleep disturbance VAS" ~ "Sleep\ndisturbance",
                                       score_lbl == "Traces of scratching" ~ "Traces\nof scratching",
                                       TRUE ~ score_lbl)
               max_score <- item_dict$Maximum[i]
               
               p <- ggplot() +
                 add_broken_pointline(poscorad, aes_x = "Day", aes_y = score_lbl, colour = "PO-SCORAD")
               if (!(score_lbl %in% subj_symp)) {
                 p <- p +
                   geom_point(data = scorad %>% rename(Score = all_of(score_lbl)),
                              aes(x = Day, y = Score, colour = "SCORAD"),
                              size = 2)
               }
               p <- p +
                 scale_colour_manual(values = palette) +
                 scale_y_continuous(limits = c(0, max_score), expand = expansion(mult = .05)) +
                 scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .02))) +
                 labs(y = score_lbl2, colour = "")

               if (score_lbl %in% intensity_signs) {
                 p <- p + theme(panel.grid.minor.y = element_blank())
               }
               
               return(p)
             })

# Remove x axis except for the last plot (unnecessary spacing)
id <- 1:(length(pl) - 1)
pl[id] <- lapply(pl[id],
                 function(x) {
                   x + theme(axis.text.x = element_blank(),
                             axis.title.x = element_blank(),
                             axis.ticks.x = element_blank())
                 })
# Legend
legend <- get_legend(pl[[1]] + theme(legend.position = "top"))
# Reorient y axis label and legend
pl <- lapply(pl,
             function(x) {
               x + theme(axis.title.y = element_text(angle = 0, vjust = 0.5),
                         legend.position = "none")
             })
# Size of each plot
sz <- item_dict %>% mutate(SizePlot = case_when(Maximum / Resolution == 100 ~ 2, TRUE ~ 1)) %>% pull(SizePlot)
sz[length(sz)]  <- sz[length(sz)] + .5

plot_grid(legend,
          plot_grid(plotlist = pl,
                    ncol = 1,
                    align = "v",
                    rel_heights = sz),
          ncol = 1, rel_heights = c(1, 8))

rhs <- plot_grid(plotlist = pl,
                 ncol = 1,
                 align = "v",
                 rel_heights = sz)

# Plot SCORAD ------------------------------------------------

p_SCORAD <- ggplot() +
  add_broken_pointline(poscorad, aes_x = "Day", aes_y = "SCORAD", colour = "PO-SCORAD") +
  geom_point(data = scorad,
             aes(x = Day, y = SCORAD, colour = "SCORAD"),
             size = 2) +
  scale_colour_manual(values = palette) +
  scale_y_continuous(limits = c(0, 103), expand = expansion(mult = c(0.02, 0))) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .02))) +
  labs(colour = "") +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  theme(legend.position = "top")

# Plot treatment ----------------------------------------------------------

treatment_names <- setNames(c("emollientCreamWithinThePast2Days", "localTreatmentWithinThePast2Days"),
                            c("Emollient\nCream", "Topical\nCorticosteroids"))

poscorad <- rename(poscorad, treatment_names)

p_treat <- lapply(names(treatment_names),
       function(y) {
         ggplot() +
           add_broken_pointline(poscorad, aes_x = "Day", aes_y = y) +
           scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .02))) +
           scale_y_continuous(breaks = c(0, 1), limits = c(0, 1)) +
           theme(axis.title.y = element_text(angle = 0, vjust = 0.5),
                 panel.grid.minor.y = element_blank())
       })

# Combine plots -----------------------------------------------------------

pl_lhs <- c(list(NULL), list(p_SCORAD), p_treat, list(NULL))

lhs <- plot_grid(plotlist = pl_lhs,
          ncol = 1, align = "v", rel_heights = c(1, 5, 1, 1, 1))

plot_grid(lhs, rhs, nrow = 1)

if (FALSE) {
  ggsave(here("results", paste0("example_data", pid, ".jpg")),
         width = 13, height = 8, units = "cm", dpi = 300, scale = 3.3, bg = "white")
}
