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
                                            TRUE ~ 0)) %>% # if alive after 5 years ~ 1, else ~ 0 %>% 
  relocate(alive_after_five_years, .after = survival_days) %>% # move column to the right of survival_days
  group_by(SEER_stage) %>%
  mutate(SEER_stage = factor(SEER_stage,
                             levels = c("In situ",
                                        "Localized",
                                        "Regional",
                                        "Distant"))) %>% # makes sorting by stage possible
  filter(tissue_type == "cancerous") # only one row per patient

# Plot: How many patients in each group (cancer stage)
text_size = 13 # for all plots
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
summarize (alive_sum = sum(alive_after_five_years), counts=n(), five_year_survival_rate = alive_sum / counts) %>%
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
  relocate(status, .after = alive_after_five_years)
surv_model <- survfit(Surv(survival_days, status) ~ SEER_stage,
                 data = data_fysr)
p_km <- ggsurvplot(surv_model,
                   data = data_fysr,
                   legend = "none",
                   font.x = text_size,
                   font.y = text_size,
                   font.tickslab = text_size)

# Putting the plots together
p_km <- p_km$plot # extract ggplot object (required for plotting together with ggplots)
box_plots <- ggarrange(p_count,
                       p_fysr,
                       labels = c("A", "B"),
                       font.label = c(size = text_size),
                       common.legend = TRUE,
                       legend = "bottom")
all_plots <- ggarrange(box_plots,
                       p_km,
                       nrow = 2,
                       heights = c(1,2),
                       labels = c("", "C"),
                       font.label = c(size = text_size))
# Add shared title
annotate_figure(all_plots,
                top = text_grob("Survival probability for different stages of cancer\n",
                                # face = "bold",
                                size = 28,
                                lineheight = 0.3))

# Save the plots as a .PNG file
ggsave("results/kaplan_meier.png",
       bg = "white",
       width = 4500,
       height = 3000,
       units = "px")

