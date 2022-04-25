library(tidyverse)
library(patchwork)

tidy_data = read_csv("./data/data_tidy.csv")

toDelete <- seq(0, length(tidy_data), 2)
data <- data[-toDelete,]

data <- data %>% 
  select(chemoradiation_therapy, survival_days, alcohol_consumption, 
         smoking, death_due_to_cancer, barretts_esophagus, nodal_involvement) %>%
  drop_na()

pl1 <- data %>%
  ggplot(aes(x = survival_days, y = smoking, colour = death_due_to_cancer)) +
  geom_point()

pl2 <- data %>%
  ggplot(aes(x = survival_days, y = alcohol_consumption, colour = death_due_to_cancer)) +
  geom_point()

pl3 <- data %>%
  ggplot(aes(x = survival_days, y = chemoradiation_therapy, colour = death_due_to_cancer)) +
  geom_point()

pl1
pl2
pl3

pl1 / pl2 / pl3



pl4 <- data %>%
  ggplot(aes(x = survival_days,
             y = smoking,
             fill = death_due_to_cancer)) +
  geom_boxplot(alpha = 0.5) +
  theme_classic(base_family = "Avenir",
                base_size = 12) +
  theme(legend.position = "none")

pl5 <- data %>%
  ggplot(aes(x = survival_days,
             y = alcohol_consumption,
             fill = death_due_to_cancer)) +
  geom_boxplot(alpha = 0.5) +
  theme_classic(base_family = "Avenir",
                base_size = 12) +
  theme(legend.position = "none")

pl6 <- data %>%
  group_by(death_due_to_cancer) %>%
  summarise(n_chemo = sum(chemoradiation_therapy == "Yes"), 
            n_smoking = sum(smoking == "Yes"),
            n_alcohol = sum(alcohol_consumption == "Yes"),
            n_barretts = sum(barretts_esophagus == "Yes"),
            n_nodal = sum(nodal_involvement == "Yes")) %>%
  pivot_longer(!death_due_to_cancer) %>%
  ggplot(aes(x = value, y = name, fill = death_due_to_cancer)) +
  geom_col(position = "dodge", alpha = 0.5) +
  theme_classic(base_family = "Avenir",
              base_size = 12) + 
  theme(legend.position = "bottom")

pl4 / pl5 / pl6
