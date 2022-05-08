library(tidyverse)
library(patchwork)
library(broom)
library(cowplot)

tidy_data_wide = read_csv("./data/data_aug_wide.csv")
tidy_data_long = read_csv("./data/data_aug_long.csv")

toDelete <- seq(0, length(tidy_data_wide), 2)
data <- tidy_data_wide[-toDelete,]

pl4 <- data %>%
  ggplot(aes(x = survival_days,
             y = smoking,
             fill = death_due_to_cancer)) +
  geom_boxplot(alpha = 0.5) +
  labs(x = "Survival days",
       y = "Smoking") +
  theme_classic(base_family = "Avenir",
                base_size = 12) +
  theme(legend.position = "none")

pl5 <- data %>%
  ggplot(aes(x = survival_days,
             y = alcohol_consumption,
             fill = death_due_to_cancer)) +
  geom_boxplot(alpha = 0.5) +
  labs(x = "Survival days",
       y = "Alcohol consumption") +
  theme_classic(base_family = "Avenir",
                base_size = 12) +
  theme(legend.position = "none")

pl6 <- data %>%
  group_by(death_due_to_cancer) %>%
  summarise(`Chemo therapy` = sum(chemoradiation_therapy == "Yes"),
            Smoking = sum(smoking == "Yes"),
            `Alcohol consumption` = sum(alcohol_consumption == "Yes"),
            `Barretts esophagus` = sum(barretts_esophagus == "Yes"),
            `Nodal involvement` = sum(nodal_involvement == "Yes")) %>%
  pivot_longer(!death_due_to_cancer) %>%
  ggplot(aes(x = value, y = name, fill = death_due_to_cancer)) +
  geom_col(position = "dodge", alpha = 0.5) +
  labs(x = "Patient count",
       y = "Conditions",
       fill = "Death due to cancer") +
  theme_classic(base_family = "Avenir",
              base_size = 12) +
  theme(legend.position = "bottom")

pl6 / (pl4 + pl5)


