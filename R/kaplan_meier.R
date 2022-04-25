library(tidyverse)
library(magrittr)
source("./R/proj_func.R")

tidy_data = read_csv("./data/data_tidy.csv")

data = tidy_data %>% 
  select(patient, tumor_type, tissue_type, cohort, barretts_esophagus,
        chemoradiation_therapy, ptnm_stage, nodal_involvement, survival_days,
        death_due_to_cancer)

data %>% 
  filter(tissue_type == "cancerous") %>% 
  mutate(survival_rates = generate_survival_probabilities(survival_days) * 100) %>%
  # group_by(chemoradiation_therapy) %>%
  mutate(prev_surv = generate_survival_probabilities(survival_days,
                                                     shift = TRUE) * 100) %>%
  pivot_longer(cols = survival_rates:prev_surv,
               names_to = "survrate_type",
               values_to = "survival_rates") %>% 
  ggplot(mapping = aes(x = survival_days,
                       y = survival_rates)) +
  geom_point() +
  geom_line()










