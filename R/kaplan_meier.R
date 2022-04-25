library(tidyverse)
library(magrittr)
source("./R/proj_func.R")

tidy_data = read_csv("./data/data_tidy.csv")

data = tidy_data %>% 
  select(patient, tumor_type, tissue_type, cohort, barretts_esophagus,
        chemoradiation_therapy, ptnm_stage, nodal_involvement, survival_days,
        death_due_to_cancer)

Kaplan_Meier_plot = data %>% 
  filter(tissue_type == "cancerous") %>% 
  mutate(survival_rates = generate_survival_probabilities(survival_days) * 100) %>%
  mutate(prev_surv = generate_survival_probabilities(survival_days,
                                                     shift = TRUE) * 100) %>%
  pivot_longer(cols = survival_rates:prev_surv,
               names_to = "survrate_type",
               values_to = "survival_rates") %>%
  mutate(survival_days = case_when(survrate_type == "prev_surv" ~ survival_days - 0.001,
                                   TRUE ~ survival_days)) %>% 
  ggplot(mapping = aes(x = survival_days,
                       y = survival_rates)) +
  xlim(0, 1825) +
  geom_line() +
  labs(title = "5-year survival",
       x = "Days",
       y = "Survival rate (%)")

ggsave(filename = "./doc/Kaplan-Meier.png",
       plot = Kaplan_Meier_plot,
       dpi = 400)
