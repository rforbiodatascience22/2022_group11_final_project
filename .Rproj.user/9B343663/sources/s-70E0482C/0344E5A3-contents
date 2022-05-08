# DEAD END: Dropped in favor of using the survival package

library("tidyverse")
library("magrittr")
source("./R/99_proj_functions.R")

tidy_data = read_csv("./data/data_aug_wide.csv")

data = tidy_data %>% 
  select(patient, tumor_type, tissue_type, cohort, barretts_esophagus,
        chemoradiation_therapy, ptnm_stage, nodal_involvement, survival_days,
        death_due_to_cancer)

data %>% 
  filter(tissue_type == "cancerous") %>% 
  mutate(survival_rates = generate_survival_probabilities(survival_days) * 100, 
         prev_surv = generate_survival_probabilities(survival_days,
                                                     shift = TRUE) * 100) %>%
  pivot_longer(cols = survival_rates:prev_surv,
               names_to = "survrate_type",
               values_to = "survival_rates") %>%
  mutate(survival_days = case_when(survrate_type == "prev_surv" ~
                                     survival_days - 0.001,
                                   TRUE ~ survival_days)) %>% 
  ggplot(mapping = aes(x = survival_days,
                       y = survival_rates)) +
  xlim(0, 1750) +
  geom_line() +
  labs(title = "5-year survival",
       x = "Days",
       y = "Survival rate (%)")

ggsave(filename = "./results/Kaplan-Meier_old.png",
       dpi = 400)
