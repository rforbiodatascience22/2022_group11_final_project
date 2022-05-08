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
  filter(tissue_type == "cancerous") # only one row per patient

# Plot: How many patients in each group (cancer stage)
text_size = 18 # for all plots
legend_name = "SEER stage"
p_count <- data_fysr %>%
  ggplot(mapping = aes(x = SEER_stage,
                       fill = SEER_stage)) +
  labs(x = "Cancer stage",
       y = "Patient count") +
  theme(legend.position = "none",
        text = element_text(size = text_size),
        axis.title.x=element_blank()) +
  scale_fill_discrete(name = legend_name) +
  geom_bar()

# Plot: 5-year survival rate per group
p_fysr <- data_fysr %>%
summarize (alive_sum = sum(alive_after_five_years),
           counts=n(),
           five_year_survival_rate = alive_sum / counts) %>%
  ggplot(mapping = aes(x = SEER_stage,
                       y = five_year_survival_rate,
                       fill = SEER_stage)) +
  theme(legend.position = "none",
        text = element_text(size = text_size),
        axis.title.x=element_blank()) +
  scale_fill_discrete(name = legend_name) +
  geom_col() +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Cancer stage",
       y = "5-year survival rate")

# Kaplan-Meier curves
data_fysr <- data_fysr %>%
  mutate(status = case_when(alive_after_five_years == 1 ~ 0,
                                       alive_after_five_years == 0 ~ 1)) %>%
  relocate(status,
           .after = alive_after_five_years)
surv_model <- survfit(Surv(survival_days,
                           status) ~ SEER_stage,
                 data = data_fysr)
p_km <- ggsurvplot(surv_model,
                   data = data_fysr,
                   legend = "none",
                   font.x = text_size,
                   font.y = text_size,
                   font.tickslab = text_size,
                   pval = TRUE)

## Putting the plots together
# Extract ggplot2 object (required for plotting together with ggplots)
p_km <- p_km$plot
box_plots <- ggarrange(p_count,
                       p_fysr,
                       labels = c("A",
                                  "B"),
                       font.label = c(size = text_size),
                       common.legend = TRUE,
                       legend = "bottom")
all_plots <- ggarrange(box_plots,
                       p_km,
                       nrow = 2,
                       heights = c(1,
                                   2),
                       labels = c("",
                                  "C"),
                       font.label = c(size = text_size))
# Add shared title
text <- "Survival probability for different stages of cancer\n"
annotate_figure(all_plots,
                top = text_grob(text,
                                size = text_size * 2.15,
                                lineheight = 0.3))

# Save the plots as a .PNG file
ggsave("results/kaplan_meier_seer.png",
       bg = "white",
       width = 4500,
       height = 3000,
       units = "px")


### Testing: Japan vs. US
surv_us_jap <- survfit(Surv(survival_days,
                            status) ~ cohort,
                       data = data_fysr)
plot_us_jap <- ggsurvplot(surv_us_jap,
                          data = data_fysr,
                          title = "All samples",
                          font.x = text_size,
                          font.y = text_size,
                          font.tickslab = text_size,
                          pval = TRUE,
                          legend.title = "Cohort",
                          legend.labs = c("Japanese",
                                          "US"))
plot_us_jap$plot <- plot_us_jap$plot +
  theme(legend.text = element_text(size = text_size),
        legend.title = element_text(size = text_size))
plot_us_jap <- plot_us_jap$plot
# US cohort contains SCC in addition to ADC. Jap does not.
data_us_jap_norm <- data_fysr %>%
  filter(tumor_type == "SCC",
         barretts_esophagus == "No")
surv_us_jap_norm <- survfit(Surv(survival_days,
                                 status) ~ cohort,
                          data = data_us_jap_norm)
plot_us_jap_norm <- ggsurvplot(surv_us_jap_norm,
                             data = data_us_jap_norm,
                             title = "SCC samples",
                             font.x = text_size,
                             font.y = text_size,
                             font.tickslab = text_size,
                             pval = TRUE)
plot_us_jap_norm <- plot_us_jap_norm$plot
us_vs_jap <- ggarrange(plot_us_jap,
                       plot_us_jap_norm,
                       common.legend = TRUE)
text_1 <- "There are no ADC samples in the Japanese cohort."
text_2 <- "Thus, the left graph is not a fair comparison."
text_all <- paste(text_1,
                  text_2,
                  sep="\n")
text_all
annotate_figure(us_vs_jap,
                top = text_grob("Survival probability by cohort\n",
                                size = text_size * 2.15,
                                lineheight = 0.3),
                bottom = text_grob(text_all,
                                   hjust = 1,
                                   x = 1,
                                   face = "italic"))
# Once filtered for SCC, the difference in outcome is dubious.
# Plot will not be included in the presentation.
ggsave("results/kaplan_meier_us_jap.png",
       bg = "white",
       width = 4500,
       height = 3000,
       units = "px")



# Testing: Which strata comparisons are statistically significant?
# data_test <- data_fysr %>%
#   filter(tumor_type == "SCC")
# surv_test <- survfit(Surv(survival_days,
#                           status) ~ smoking,
#                      data = data_test)
# plot_test <- ggsurvplot(surv_test,
#                         data = data_test,
#                         pval = TRUE)
# plot_test

# p-values for different comparisons (all samples):
# death_due_to_cancer < 0.0001
# barretts_esophagus = 0.0099
# SEER_stage = 0.021
# ptnm_stage = 0.025
# tumor_type = 0.037
# ------------------------------------------------------------------------------
# nodal_involvement = 0.073
# chemoradiation_therapy = 0.075
# alcohol_consumption = 0.38
# smoking = 0.98
# 
# ______________________________________________________________________________
# # p-values for different comparisons (ADC only):
# death_due_to_cancer < 0.0001
# ------------------------------------------------------------------------------
# ptnm_stage = 0.49
# SEER_stage = 0.73
# barretts_esophagus = 0.074
# chemoradiation_therapy = 0.11
# nodal_involvement = 0.32
# alcohol_consumption = 0.53
# smoking = 0.97
# 
# ______________________________________________________________________________
# # p-values for different comparisons (SCC only):
# death_due_to_cancer = 0.0075
# SEER_stage = 0.016
# ptnm_stage = 0.027
# ------------------------------------------------------------------------------
# chemoradiation_therapy = 0.12
# nodal_involvement = 0.26
# alcohol_consumption = 0.52
# smoking = 0.67

## Plot the four significant strata comparisons
# Barrett's esophagus
surv_barr <- survfit(Surv(survival_days,
                          status) ~ barretts_esophagus,
                     data = data_fysr)
plot_barr <- ggsurvplot(surv_barr,
                        data = data_fysr,
                        font.x = text_size,
                        font.y = text_size,
                        font.tickslab = text_size,
                        pval = TRUE,
                        legend.labs = c("Absent",
                                        "Present"),
                        legend.title = "Barrett's esophagus"
                        )
plot_barr <- plot_barr$plot +
  theme(legend.text = element_text(size = text_size),
        legend.title = element_text(size = text_size))

# SEER stage
surv_seer <- survfit(Surv(survival_days,
                          status) ~ SEER_stage,
                     data = data_fysr)
plot_seer <- ggsurvplot(surv_seer,
                        data = data_fysr,
                        font.x = text_size,
                        font.y = text_size,
                        font.tickslab = text_size,
                        pval = TRUE,
                        legend.labs = c("In situ",
                        "Localized",
                        "Regional",
                        "Distant"),
                        legend.title = "SEER stage"
)
plot_seer <- plot_seer$plot +
  theme(legend.text = element_text(size = text_size),
        legend.title = element_text(size = text_size))

# pTNM stage
surv_ptnm <- survfit(Surv(survival_days,
                          status) ~ ptnm_stage,
                     data = data_fysr)
plot_ptnm <- ggsurvplot(surv_ptnm,
                        data = data_fysr,
                        font.x = text_size,
                        font.y = text_size,
                        font.tickslab = text_size,
                        pval = TRUE,
                        legend.labs = c("0",
                                        "I",
                                        "II",
                                        "IIA",
                                        "IIB",
                                        "III",
                                        "IV"),
                        legend.title = "pTNM stage"
)
plot_ptnm <- plot_ptnm$plot +
  theme(legend.text = element_text(size = text_size),
        legend.title = element_text(size = text_size))

# Tumor type
surv_tumo <- survfit(Surv(survival_days,
                          status) ~ tumor_type,
                     data = data_fysr)
plot_tumo <- ggsurvplot(surv_tumo,
                        data = data_fysr,
                        font.x = text_size,
                        font.y = text_size,
                        font.tickslab = text_size,
                        pval = TRUE,
                        legend.title = "Tumor type",
                        legend.labs = c("ADC",
                                        "SCC")
)
plot_tumo <- plot_tumo$plot +
  theme(legend.text = element_text(size = text_size),
        legend.title = element_text(size = text_size))

# Plot all of them together
all_strata <- ggarrange(plot_tumo,
                        plot_barr,
                        plot_seer,
                        plot_ptnm)
text <- "Visualization of all significantly different strata\n"
annotate_figure(all_strata,
                top = text_grob(text,
                                size = text_size * 2.15,
                                lineheight = 0.3))
ggsave("results/kaplan_meier_all.png",
       bg = "white",
       width = 4500,
       height = 3000,
       units = "px")

## Decision: Make a combined plot of SEER stage and tumor type
## Include a p_count plot for each strata (how many in each group)
# P_count plot for tumor type
p_count_tumo <- data_fysr %>%
  ggplot(mapping = aes(x = tumor_type,
                       fill = tumor_type)) +
  labs(x = "Tumor type",
       y = "Patient count") +
  theme(legend.position = "bottom",
        text = element_text(size = text_size),
        axis.title.x=element_blank()) +
  scale_fill_discrete(name = "Tumor type") +
  geom_bar()
# Kaplan-Meier plot for the tumor type
p_tumo <- ggsurvplot(surv_tumo,
                     data = data_fysr,
                     font.x = text_size,
                     font.y = text_size,
                     font.tickslab = text_size,
                     pval = TRUE,
                     legend = "none"
)
p_tumo <- p_tumo$plot +
  theme(legend.text = element_text(size = text_size),
        legend.title = element_text(size = text_size))

# Combined tumor type plot
tumo_plots <- ggarrange(p_count_tumo,
                        p_tumo,
                        nrow = 2,
                        heights = c(1,
                                    2),
                        font.label = c(size = text_size))

# P_count and Kaplan-Meier plot for the SEER stage
p_count_seer <- p_count +
  theme(legend.position = "bottom")

# Combined SEER stage plot
seer_plots <- ggarrange(p_count_seer,
                        p_km,
                        nrow = 2,
                        heights = c(1,
                                    2),
                        font.label = c(size = text_size))

# Combine all plots
tumo_seer_plots <- ggarrange(tumo_plots,
                             seer_plots)

# Add shared title
annotate_figure(tumo_seer_plots,
                top = text_grob("Survival probability for different groups\n",
                                size = text_size * 2.15,
                                lineheight = 0.3))

# Save the plots as a .PNG file
ggsave("results/kaplan_meier_tumo_seer.png",
       bg = "white",
       width = 5000,
       height = 2500,
       units = "px")
