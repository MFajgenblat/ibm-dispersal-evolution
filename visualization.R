#-------------------------------------------------------------------------------
# Loading packages
#-------------------------------------------------------------------------------

library(tidyverse)
library(ggh4x)
library(patchwork)
library(ggtext)

#-------------------------------------------------------------------------------
# Reading simulation results
#-------------------------------------------------------------------------------

results_general_scenarios <- read.csv("results_general_scenarios.csv", sep=";")
results_contrasting_scenarios <- read.csv("results_contrasting_scenarios.csv", sep=";")

#-------------------------------------------------------------------------------
# Figure 1: Dispersal phenotypes across treatments
#-------------------------------------------------------------------------------

results_general_scenarios %>%
  filter(evodisp_scenario == 1,
         delta_fix == 0.05,
         Generation > 1900) %>%
  mutate(species_number = case_when(species_pool == 1 ~ 1,
                                    species_pool == 2 ~ 2,
                                    species_pool == 3 ~ 5),
         lambda = factor(case_when(lambda == 0.1 ~ "Low",
                                   lambda == 1 ~ "Moderate",
                                   lambda == 10 ~ "High"),
                         levels = c("Low", "Moderate", "High")),
         chi = factor(case_when(chi == 0 ~ "No cost of dispersal",
                                chi == 1 ~ "Cost of dispersal"),
                      levels = c("No cost of dispersal", "Cost of dispersal")),
         psi = factor(case_when(psi == 0 ~ "No patch\nextinction",
                                psi == 0.01 ~ "Periodic patch\nextinction"))) %>%
  group_by(species_number, lambda, chi, psi, replicate_id) %>%
  summarise(Dispersal_phenotype = mean(Dispersal_phenotype)) %>%
  ggplot(aes(x = factor(species_number), y = Dispersal_phenotype/5, color = factor(lambda), fill = factor(lambda))) +
  stat_summary(position = position_dodge(width = 0.75),fun.data = "mean_sdl", fun.args = list(mult = 1), size = 0.06, linewidth = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), shape = 16, size = 0.5, alpha = 0.25) +
  scale_x_discrete("Regional species richness") +
  scale_y_continuous("Dispersal propensity") +
  scale_color_manual("Mainland immigration", values = c("#3b6f75", "#a32f66", "#8c7837"),
                     guide = guide_legend(override.aes = list(linetype = "blank", size = 0.4))) +
  scale_fill_manual("Mainland immigration", values = c("#3b6f75", "#a32f66", "#8c7837")) +
  facet_nested(. ~ factor(chi) + factor(psi), nest_line = element_line(linetype = 1)) +
  theme(legend.position = "bottom",
        legend.key = element_blank(),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7),
        legend.margin = margin(0,0,0,0),
        legend.key.size = unit(0.4, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8),
        ggh4x.facet.nestline = element_line(size = 0.15),
        axis.title = element_text(face = "bold", size = 8),
        axis.text = element_text(size = 7),
        panel.spacing.y = unit(0.2, "cm"),
        panel.background = element_blank(),
        panel.grid = element_line(color = "grey93"),
        panel.border = element_rect(fill = NA, color = "black"))
ggsave("Figure_1.png", width = 16, height = 7, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Figure 2: Monopolization index
#-------------------------------------------------------------------------------

MI <- results_general_scenarios %>%
  filter(Generation > 1000) %>%
  group_by(lambda, chi, psi, replicate_id, evodisp_scenario, species_pool, delta_fix) %>%
  summarise(MI = 1.25*mean((Patch != Species) * (1 - abs(environment[Patch] - Environmental_phenotype)))) %>%
  mutate(species_number = case_when(species_pool == 1 ~ 1,
                                    species_pool == 2 ~ 2,
                                    species_pool == 3 ~ 5),
         lambda = factor(case_when(lambda == 0.1 ~ "Low seeding rate from mainland",
                                   lambda == 1 ~ "Moderate seeding rate from mainland",
                                   lambda == 10 ~ "High seeding rate from mainland"),
                         levels = c("Low seeding rate from mainland", "Moderate seeding rate from mainland", "High seeding rate from mainland")),
         chi = factor(case_when(chi == 0 ~ "No cost of dispersal",
                                chi == 1 ~ "Cost of dispersal"),
                      levels = c("No cost of dispersal", "Cost of dispersal")),
         psi = factor(case_when(psi == 0 ~ "No patch\nextinction",
                                psi == 0.01 ~ "Periodic patch\nextinction"),
                      levels = c("No patch\nextinction", "Periodic patch\nextinction")),
         species_number = paste(species_number, "species"),
         dispersal = factor(case_when(evodisp_scenario == 1 ~ "Evo",
                                      evodisp_scenario == 2 ~ as.character(delta_fix/5)),
                            levels = c("0.01", "Evo", "0.25")),
         treatment = paste(species_number, lambda, chi, psi, sep = "\n")) %>%
  ungroup() %>%
  mutate(treatment = fct_reorder(treatment, MI, mean))
MI_A <- MI %>%
  ggplot(aes(x = treatment, y = MI, color = dispersal)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.5), shape = 16, size = 0.4, alpha = 0.5) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), size = 0.005, linewidth = 0.15, position = position_dodge(width = 0.5)) +
  scale_y_continuous("Monopolization index", breaks = seq(0, 1, by = 0.2)) +
  scale_color_manual("Dispersal propensity", values = c("#075C73", "#E29507", "#CD2930"),
                     guide = guide_legend(override.aes = list(linetype = "blank", size = 0.4))) +
  coord_cartesian(clip = "off", ylim = c(0,1)) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        panel.background = element_blank(),
        panel.grid = element_line(color = "grey93"),
        panel.border = element_rect(fill = NA, color = "black", size = 0.2),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 8),
        axis.text.y = element_text(size = 7),
        legend.position = "top",
        legend.key = element_blank(),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7))
MI_B <- MI %>%
  select(species_number, lambda, chi, psi, treatment) %>%
  distinct() %>%
  mutate(species_number = as.character(case_when(species_number == "1 species" ~ 1,
                                                 species_number == "2 species" ~ 2,
                                                 species_number == "5 species" ~ 5)),
         lambda = case_when(lambda == "Low seeding rate from mainland" ~ "Low",
                            lambda == "Moderate seeding rate from mainland" ~ "Mod.",
                            lambda == "High seeding rate from mainland" ~ "High"),
         chi = case_when(chi == "No cost of dispersal" ~ "No",
                         chi == "Cost of dispersal" ~ "Yes"),
         psi = case_when(psi == "No patch\nextinction" ~ "No",
                         psi == "Periodic patch\nextinction" ~ "Yes")) %>%
  pivot_longer(c(species_number, lambda, chi, psi)) %>%
  mutate(name = case_when(name == "species_number" ~ "Number of\nspecies",
                          name == "lambda" ~ "Mainland\nimmigration",
                          name == "chi" ~ "Cost of\ndispersal",
                          name == "psi" ~ "Periodic patch\nextinction")) %>%
  ggplot() +
  geom_text(aes(x = treatment, y = name, label = value), size = 1.4, color = "grey40") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 7))
MI_A + MI_B + plot_layout(ncol = 1, heights = c(2,1))
ggsave("Figure_2.pdf", width = 16, height = 12, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Figure 3: Contrasting scenarios - population density
#-------------------------------------------------------------------------------

results_contrasting_scenarios %>%
  filter(Generation >= 1000) %>%
  group_by(replicate_id, lambda, psi, evodisp_scenario, delta_fix, Species) %>%
  summarize(MI = mean(Dispersal_phenotype)/5) %>%
  mutate(psi = factor(case_when(psi == 0 ~ "No patch\nextinction",
                                psi == 0.01 ~ "Periodic patch\nextinction"),
                      levels = c("No patch\nextinction", "Periodic patch\nextinction")),
         lambda = factor(case_when(lambda == 0.1 ~ "Low mainland dispersal",
                                   lambda == 1 ~ "Moderate mainland dispersal",
                                   lambda == 10 ~ "High mainland dispersal"),
                         levels = c("Low mainland dispersal", "Moderate mainland dispersal", "High mainland dispersal")),
         MI = case_when(evodisp_scenario == 2 ~ delta_fix/5,
                        (evodisp_scenario == 3 & Species %in% c(1)) ~ delta_fix/5,
                        (evodisp_scenario == 4 & Species %in% c(3)) ~ delta_fix/5,
                        (evodisp_scenario == 5 & Species %in% c(1,2,4,5)) ~ delta_fix/5,
                        (evodisp_scenario == 6 & Species %in% c(2,3,4,5)) ~ delta_fix/5,
                        T ~ MI),
         evodisp_scenario = factor(c("All", "None", "All but species 1", "All but species 3", "Only species 3", "Only species 1")[evodisp_scenario],
                                   levels = c("All", "All but species 1", "All but species 3", "Only species 3", "Only species 1", "None")),
         delta_fix = factor(case_when(delta_fix == 0.05 ~ "Low fixed dispersal propensity",
                                      delta_fix == 1.25 ~ "High fixed dispersal propensity"),
                            levels = c("Low fixed dispersal propensity", "High fixed dispersal propensity"))) %>%
  filter(!c(evodisp_scenario %in% c("All but species 1", "All but species 3"))) %>%
  ggplot(aes(x = evodisp_scenario, y = MI, color = factor(Species))) +
  stat_summary(position = position_dodge(width = 0.75),fun.data = "mean_sdl", fun.args = list(mult = 1), size = 0.3, stroke = NA, linewidth = 0.3) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), shape = 16, size = 0.2, alpha = 0.5) +
  scale_x_discrete("Species with an evolvable dispersal propensity") +
  scale_y_continuous("Average dispersal phenotype") +
  scale_color_brewer("Species", palette = "Dark2",
                     guide = guide_legend(override.aes = list(linetype = "blank", size = 0.65))) +
  facet_nested(delta_fix + psi ~ lambda, nest_line = element_line(linetype = 1)) +
  theme(legend.position = "bottom",
        legend.key = element_blank(),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7),
        legend.margin = margin(0,0,0,0),
        legend.key.size = unit(0.4, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8),
        strip.placement = "outside",
        ggh4x.facet.nestline = element_line(size = 0.15),
        axis.title = element_text(face = "bold", size = 8),
        axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 7),
        panel.spacing.y = unit(0.2, "cm"),
        panel.background = element_blank(),
        panel.grid = element_line(color = "grey93"),
        panel.border = element_rect(fill = NA, color = "black"))
ggsave("Figure_3.pdf", width = 16, height = 16, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Figure 4: Contrasting scenarios - monopolization index
#-------------------------------------------------------------------------------

results_contrasting_scenarios %>%
  filter(Generation >= 1000) %>%
  group_by(replicate_id, lambda, psi, evodisp_scenario, delta_fix, Species) %>%
  summarize(MI = 1.25*mean((Patch != Species) * (1 - abs(environment[Patch] - Environmental_phenotype)))) %>%
  mutate(psi = factor(case_when(psi == 0 ~ "No patch\nextinction",
                                psi == 0.01 ~ "Periodic patch\nextinction"),
                      levels = c("No patch\nextinction", "Periodic patch\nextinction")),
         lambda = factor(case_when(lambda == 0.1 ~ "Low mainland dispersal",
                                   lambda == 1 ~ "Moderate mainland dispersal",
                                   lambda == 10 ~ "High mainland dispersal"),
                         levels = c("Low mainland dispersal", "Moderate mainland dispersal", "High mainland dispersal")),
         evodisp_scenario = factor(c("All", "None", "All but species 1", "All but species 3", "Only species 3", "Only species 1")[evodisp_scenario],
                                   levels = c("All", "All but species 1", "All but species 3", "Only species 3", "Only species 1", "None")),
         delta_fix = factor(case_when(delta_fix == 0.05 ~ "Low fixed dispersal propensity",
                                           delta_fix == 1.25 ~ "High fixed dispersal propensity"),
                                 levels = c("Low fixed dispersal propensity", "High fixed dispersal propensity"))) %>%
  filter(!c(evodisp_scenario %in% c("All but species 1", "All but species 3"))) %>%
  ggplot(aes(x = evodisp_scenario, y = MI, color = factor(Species))) +
  stat_summary(position = position_dodge(width = 0.75),fun.data = "mean_sdl", fun.args = list(mult = 1), size = 0.3, stroke = NA, linewidth = 0.3) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), shape = 16, size = 0.2, alpha = 0.5) +
  scale_x_discrete("Species with an evolvable dispersal propensity") +
  scale_y_continuous("Monopolization index", breaks = seq(0, 1, by = 0.2)) +
  scale_color_brewer("Species", palette = "Dark2",
                     guide = guide_legend(override.aes = list(linetype = "blank", size = 0.65))) +
  coord_cartesian(ylim = c(0,1)) +
  facet_nested(delta_fix + psi ~ lambda, nest_line = element_line(linetype = 1)) +
  theme(legend.position = "bottom",
        legend.key = element_blank(),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7),
        legend.margin = margin(0,0,0,0),
        legend.key.size = unit(0.4, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8),
        strip.placement = "outside",
        ggh4x.facet.nestline = element_line(size = 0.15),
        axis.title = element_text(face = "bold", size = 8),
        axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 7),
        panel.spacing.y = unit(0.2, "cm"),
        panel.background = element_blank(),
        panel.grid = element_line(color = "grey93"),
        panel.border = element_rect(fill = NA, color = "black", size = 0.4))
ggsave("Figure_4.pdf", width = 16, height = 16, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Figure S1: Contrasting scenarios - population density
#-------------------------------------------------------------------------------

results_contrasting_scenarios %>%
  filter(Generation >= 1000) %>%
  group_by(replicate_id, lambda, psi, evodisp_scenario, delta_fix, Species) %>%
  summarize(MI = n()/10) %>%
  mutate(psi = factor(case_when(psi == 0 ~ "No patch\nextinction",
                                psi == 0.01 ~ "Periodic patch\nextinction"),
                      levels = c("No patch\nextinction", "Periodic patch\nextinction")),
         lambda = factor(case_when(lambda == 0.1 ~ "Low mainland dispersal",
                                   lambda == 1 ~ "Moderate mainland dispersal",
                                   lambda == 10 ~ "High mainland dispersal"),
                         levels = c("Low mainland dispersal", "Moderate mainland dispersal", "High mainland dispersal")),
         evodisp_scenario = factor(c("All", "None", "All but species 1", "All but species 3", "Only species 3", "Only species 1")[evodisp_scenario],
                                   levels = c("All", "All but species 1", "All but species 3", "Only species 3", "Only species 1", "None")),
         delta_fix = factor(case_when(delta_fix == 0.05 ~ "Low fixed dispersal propensity",
                                           delta_fix == 1.25 ~ "High fixed dispersal propensity"),
                                 levels = c("Low fixed dispersal propensity", "High fixed dispersal propensity"))) %>%
  filter(!c(evodisp_scenario %in% c("All but species 1", "All but species 3"))) %>%
  ggplot(aes(x = evodisp_scenario, y = MI, color = factor(Species))) +
  stat_summary(position = position_dodge(width = 0.75),fun.data = "mean_sdl", fun.args = list(mult = 1), size = 0.3, stroke = NA, linewidth = 0.3) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), shape = 16, size = 0.2, alpha = 0.5) +
  scale_x_discrete("Species with an evolvable dispersal propensity") +
  scale_y_continuous("Average population size per generation") +
  scale_color_brewer("Species", palette = "Dark2",
                     guide = guide_legend(override.aes = list(linetype = "blank", size = 0.65))) +
  facet_nested(delta_fix + psi ~ lambda, nest_line = element_line(linetype = 1)) +
  theme(legend.position = "bottom",
        legend.key = element_blank(),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7),
        legend.margin = margin(0,0,0,0),
        legend.key.size = unit(0.4, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8),
        strip.placement = "outside",
        ggh4x.facet.nestline = element_line(size = 0.15),
        axis.title = element_text(face = "bold", size = 8),
        axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 7),
        panel.spacing.y = unit(0.2, "cm"),
        panel.background = element_blank(),
        panel.grid = element_line(color = "grey93"),
        panel.border = element_rect(fill = NA, color = "black"))
ggsave("Figure_S1.pdf", width = 16, height = 16, units = "cm", dpi = 600)

#-------------------------------------------------------------------------------
# Figure S2: Plotting individual replicates for each scenario
#-------------------------------------------------------------------------------

plotting_grid_original <- expand.grid(species_pool = 1:3,
                                      evodisp_scenario = 1,
                                      delta_fix = 0.05,
                                      lambda = c(0.1,1,10),
                                      chi = c(0,1),
                                      psi = c(0,0.01),
                                      replicate_id = 1)
plotting_grid_renamed <- plotting_grid_original %>%
  mutate(species_number = case_when(species_pool == 1 ~ 1,
                                    species_pool == 2 ~ 2,
                                    species_pool == 3 ~ 5),
         species_number = paste(species_number, "species"),
         lambda = case_when(lambda == 0.1 ~ "low seeding rate from mainland",
                            lambda == 1 ~ "moderate seeding rate from mainland",
                            lambda == 10 ~ "high seeding rate from mainland"),
         chi = case_when(chi == 0 ~ "no cost of dispersal",
                         chi == 1 ~ "cost of dispersal"),
         psi = case_when(psi == 0 ~ "no patch extinction",
                         psi == 0.01 ~ "periodic patch extinction"))
for (i in 1:nrow(plotting_grid_original)) {
  results_general_scenarios %>%
    filter(evodisp_scenario == 1) %>%
    filter(species_pool == plotting_grid_original$species_pool[i], lambda == plotting_grid_original$lambda[i], chi == plotting_grid_original$chi[i], psi == plotting_grid_original$psi[i], replicate_id == 1) %>%
    mutate(species_number = case_when(species_pool == 1 ~ 1,
                                      species_pool == 2 ~ 2,
                                      species_pool == 3 ~ 5),
           species_number = paste(species_number, "species"),
           lambda = case_when(lambda == 0.1 ~ "low seeding rate from mainland",
                              lambda == 1 ~ "moderate seeding rate from mainland",
                              lambda == 10 ~ "high seeding rate from mainland"),
           chi = case_when(chi == 0 ~ "no cost of dispersal",
                           chi == 1 ~ "cost of dispersal"),
           psi = case_when(psi == 0 ~ "no patch extinction",
                           psi == 0.01 ~ "periodic patch extinction")) %>%
    sample_n(10000) %>%
    pivot_longer(c(Environmental_phenotype, Dispersal_phenotype)) %>%
    ggplot() +
    geom_point(aes(x = Generation, y = value, color = factor(Species)), alpha = 0.25, shape = 16, size = 0.2) +
    scale_color_brewer("Species", palette = "Dark2", guide = "none") +
    scale_y_continuous("Phenotype", limits = c(0,1), expand = c(0,0), breaks = seq(0, 1, by = 0.2)) +
    facet_grid(gsub("_phenotype", "", name) ~ paste("Patch", Patch), switch = "y") +
    ggtitle(paste0(c(letters, apply(expand.grid(letters, letters)[,c(2,1)], 1, paste0, collapse = ""))[i], ") ", plotting_grid_renamed$species_number[i], ", ", plotting_grid_renamed$lambda[i], ", ", plotting_grid_renamed$chi[i], ", ", plotting_grid_renamed$psi[i])) +
    theme(plot.title = element_text(face = "bold", size = 8),
          legend.position = "bottom",
          legend.key = element_blank(),
          legend.title = element_text(face = "bold", size = 8),
          legend.text = element_text(size = 7),
          legend.margin = margin(0,0,0,0),
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 8),
          axis.title = element_text(face = "bold", size = 8),
          axis.text = element_text(size = 7),
          panel.spacing.y = unit(0.5, "cm"),
          panel.background = element_blank(),
          panel.grid = element_line(color = "grey93"),
          panel.border = element_rect(fill = NA, color = "black"))
  ggsave(paste0("Figure_S2", c(letters, apply(expand.grid(letters, letters)[,c(2,1)], 1, paste0, collapse = ""))[i], ".png"), width = 16, height = 7, units = "cm", dpi = 600)
}
