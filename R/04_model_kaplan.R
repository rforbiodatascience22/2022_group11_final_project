library(tidyverse)
library(gridExtra)
library(survival)
library(survminer)

# Read the data
data <- read_csv("data/data_aug_wide.csv",
                 show_col_types = FALSE) %>%
  as_tibble()

# Make the data suitable for plotting 5-year survival rates
data_fysr <- data %>%
  mutate(alive_after_five_years = case_when(survival_days == "1826" ~ 1,
                                            TRUE ~ 0)) %>%
  relocate(alive_after_five_years,
           .after = survival_days) %>%
  group_by(SEER_stage) %>%
  # Make sorting by SEER stage possible:
  mutate(SEER_stage = factor(SEER_stage,
                             levels = c("In situ",
                                        "Localized",
                                        "Regional",
                                        "Distant"))) %>%
  filter(tissue_type == "cancerous") %>%  # only one row per patient
  # add alive status (required for survival package)
  mutate(status = case_when(alive_after_five_years == 1 ~ 0,
                            alive_after_five_years == 0 ~ 1)) %>%
  relocate(status,
           .after = alive_after_five_years)

# Set the text size for all plots
text_size = 18

# Plot the group counts for each tumor type
count_tumor <- data_fysr %>%
  ggplot(mapping = aes(x = tumor_type,
                       fill = tumor_type)) +
  labs(x = "Tumor type",
       y = "Patient count") +
  theme(legend.position = "bottom",
        text = element_text(size = text_size),
        axis.title.x=element_blank()) +
  scale_fill_discrete(name = "Tumor type") +
  geom_bar()

# Plot the Kaplan-Meier curves for each tumor type
surv_tumor <- survfit(Surv(survival_days,
                          status) ~ tumor_type,
                     data = data_fysr)
kaplan_tumor <- ggsurvplot(surv_tumor,
                     data = data_fysr,
                     font.x = text_size,
                     font.y = text_size,
                     font.tickslab = text_size,
                     pval = TRUE,
                     legend = "none")
kaplan_tumor <- kaplan_tumor$plot
kaplan_tumor <- kaplan_tumor +
  theme(legend.text = element_text(size = text_size),
        legend.title = element_text(size = text_size))

# Combine the tumor type plots
tumor_plots <- ggarrange(count_tumor,
                        kaplan_tumor,
                        nrow = 2,
                        heights = c(1,
                                    2),
                        font.label = c(size = text_size))

# Plot the group counts for each SEER stage
count_seer <- data_fysr %>%
  ggplot(mapping = aes(x = SEER_stage,
                       fill = SEER_stage)) +
  labs(x = "SEER stage",
       y = "Patient count") +
  theme(legend.position = "bottom",
        text = element_text(size = text_size),
        axis.title.x=element_blank()) +
  scale_fill_discrete(name = "SEER stage") +
  geom_bar()

# Plot the Kaplan-Meier curves for each SEER stage
surv_seer <- survfit(Surv(survival_days,
                           status) ~ SEER_stage,
                      data = data_fysr)
kaplan_seer <- ggsurvplot(surv_seer,
                           data = data_fysr,
                           font.x = text_size,
                           font.y = text_size,
                           font.tickslab = text_size,
                           pval = TRUE,
                           legend = "none")
kaplan_seer <- kaplan_seer$plot
kaplan_seer <- kaplan_seer +
  theme(legend.text = element_text(size = text_size),
        legend.title = element_text(size = text_size))

# Combine the SEER stage plots
seer_plots <- ggarrange(count_seer,
                        kaplan_seer,
                        nrow = 2,
                        heights = c(1,
                                    2),
                        font.label = c(size = text_size))

# Combine all plots
tumor_seer_plots <- ggarrange(tumor_plots,
                             seer_plots)

# Add shared title
annotate_figure(tumor_seer_plots,
                top = text_grob("Survival probability for different groups\n",
                                size = text_size * 2.15,
                                lineheight = 0.3))

# Save the plots as a .PNG file
ggsave("results/kaplan_meier_tumo_seer.png",
       bg = "white",
       width = 5000,
       height = 2500,
       units = "px")