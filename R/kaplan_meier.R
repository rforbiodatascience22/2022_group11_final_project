library(tidyverse)
library(magrittr)
source("./R/proj_func.R")

tidy_data = read_csv("./data/data_tidy.csv")

data = tidy_data %>% 
  select(patient, tumor_type, tissue_type, cohort, barretts_esophagus,
        chemoradiation_therapy, ptnm_stage, nodal_involvement, survival_days,
        death_due_to_cancer)

data %>% 
  filter(tissue_type == "cancerous",
         death_due_to_cancer == "Yes") %>% 
  mutate(survival_rates = generate_survival_probabilities(survival_days) * 100) %>%
  group_by(chemoradiation_therapy)
  ggplot(mapping = aes(x = survival_days,
                       y = survival_rates,
                       color = chemoradiation_therapy)) +
  geom_point() +
  geom_line()






