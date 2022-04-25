library(tidyverse)

data_tidy <- read_csv("data/data_tidy.csv") %>% 
  as_tibble()

# how many in each group
data_tidy %>%
  ggplot(mapping = aes(x = SEER_stage)) +
  geom_bar()

# plot
data_tidy %>% 
  mutate(five_year_survival_rate = case_when(survival_days == "1826" ~ 1,
                                          TRUE ~ 0)) %>% # if alive after 5 years ~ 1
  relocate(five_year_survival_rate, .after = survival_days) %>% 
  group_by(SEER_stage) %>% 
  summarize (alive = sum(five_year_survival_rate), counts=n(), fysr = alive / counts) %>%
  mutate(SEER_stage = factor(SEER_stage,
                             levels = c("In situ",
                                        "Localized",
                                        "Regional",
                                        "Distant"))) %>%
  ggplot(mapping = aes(x = SEER_stage,
                       y = fysr)) +
  geom_col()
#  View()
