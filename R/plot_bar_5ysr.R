library(tidyverse)
library(gridExtra)
library(survival)
library(survminer)

#NOTE: "./data/data_tidy.csv" IS NOW "./data/data_aug_wide.csv".
#YOU SHOULD CHANGE IT IN YOUR READ_CSV FUNCTION
data_tidy <- read_csv("data/data_tidy.csv") %>% 
  as_tibble()

# Calculating the 5-year survival rate
data_fysr <- data_tidy %>% 
  mutate(alive_after_five_years = case_when(survival_days == "1826" ~ 1,
                                            TRUE ~ 0)) %>% # if alive after 5 years ~ 1 %>% 
  relocate(alive_after_five_years, .after = survival_days) %>% # move column to the right of survival_days
  group_by(SEER_stage) %>%
  mutate(SEER_stage = factor(SEER_stage,
                             levels = c("In situ",
                                        "Localized",
                                        "Regional",
                                        "Distant"))) %>% # makes sorting by stage possible
  filter(tissue_type == "cancerous") # only one row per patient

# Plot: How many patients in each group (cancer stage)
p_count <- data_fysr %>%
  ggplot(mapping = aes(x = SEER_stage,
                       fill = SEER_stage)) +
  labs(x = "Cancer stage",
       y = "Patient count") +
  theme(legend.position = "none") +
  geom_bar()

# Plot: 5-year survival rate per group
p_fysr <- data_fysr %>%
summarize (alive_sum = sum(alive_after_five_years), counts=n(), five_year_survival_rate = alive_sum / counts) %>%
  ggplot(mapping = aes(x = SEER_stage,
                       y = five_year_survival_rate,
                       fill = SEER_stage)) +
  theme(legend.position = "none") +
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
                   legend.labs = c("In situ", "Localized", "Regional", "Distant"),
                   legend.title = "SEER stage")
### MEMO TO SELF: Include the legend from the bar plot instead of the ggsurv
### Reason: The legend symbols will be boxes, which generalizes better to all of the plots
# Putting the plots together
grid.arrange(p_count, p_fysr, nrow = 1)
p_km <- p_km$plot # extract ggplot object (required for plotting together with ggplots)

gl = list()
gl[[1]] = p_count
gl[[2]] = p_fysr
gl[[3]] = p_km
grid.arrange(grobs = gl,
             heights = c(1,2),
             layout_matrix = rbind(c(1,2),
                                   rbind(c(3,3))))