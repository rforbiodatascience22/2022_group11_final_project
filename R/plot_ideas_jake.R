library(tidyverse)
library(patchwork)

#NOTE: "./data/data_tidy.csv" IS NOW "./data/data_aug_wide.csv".
#YOU SHOULD CHANGE IT IN YOUR READ_CSV FUNCTION
tidy_data_wide = read_csv("./data/data_aug_wide.csv")

toDelete <- seq(0, length(tidy_data_wide), 2)
data <- tidy_data_wide[-toDelete,]

data <- data %>% 
  select(chemoradiation_therapy, survival_days, alcohol_consumption, 
         smoking, death_due_to_cancer, barretts_esophagus, nodal_involvement) %>%
  drop_na()

pl1 <- data %>%
  ggplot(aes(x = survival_days, y = smoking, fill = death_due_to_cancer)) +
  geom_boxplot(alpha = 0.5) +
  theme_classic(base_family = "Avenir",
                base_size = 12) +
  theme(legend.position = "none")

pl2 <- data %>%
  ggplot(aes(x = survival_days, y = alcohol_consumption, fill = death_due_to_cancer)) +
  geom_boxplot(alpha = 0.5) +
  theme_classic(base_family = "Avenir",
                base_size = 12) +
  theme(legend.position = "none")

pl3 <- data %>%
  ggplot(aes(x = survival_days, y = chemoradiation_therapy, fill = death_due_to_cancer)) +
  geom_boxplot(alpha = 0.5) +
  theme_classic(base_family = "Avenir",
                base_size = 12) +
  theme(legend.position = "bottom")

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

# PCA
tidy_data_long = read_csv("./data/data_aug_long.csv")

library(broom)
library(cowplot)

pca_fit <- tidy_data_long %>% 
  select(patient, Expression, Mature_MiRNA) %>%
  drop_na() %>%
  group_by(Mature_MiRNA, patient) %>%
  mutate(Mean = mean(Expression)) %>%
  select(-Expression) %>%
  distinct() %>%
  pivot_wider(names_from = Mature_MiRNA, values_from = Mean) %>%
  select_if(~ !any(is.na(.)))%>%
  ungroup() %>%
  select(-patient) %>%
  scale() %>%
  prcomp(scale = TRUE)

pl7 <- pca_fit %>%
  augment(tidy_data_wide) %>% # add original dataset back in
  ggplot(aes(.fittedPC1, .fittedPC2, color = SEER_stage)) + 
  geom_point(size = 1.5) +
  scale_color_manual(
    values = c(Distant = "#00b300", Localized = "#0072B2", Regional = "#FF0000",
               `In situ` = "#800080")
  ) +
  theme_half_open(12) + background_grid()

# define arrow style for plotting
arrow_style <- arrow(
  angle = 20, ends = "first", type = "closed", length = grid::unit(8, "pt")
)

# plot rotation matrix
pl8 <- pca_fit %>%
  tidy(matrix = "rotation") %>%
  pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value") %>%
  ggplot(aes(PC1, PC2)) +
  geom_segment(xend = 0, yend = 0, arrow = arrow_style) +
  geom_text(
    aes(label = column),
    hjust = 1, nudge_x = -0.02, 
    color = "#904C2F", size = 2,
    check_overlap = T
  ) +
  xlim(-0.5, .5) + ylim(-.5, 0.5) +
  coord_fixed() + # fix aspect ratio to 1:1
  theme_minimal_grid(12)

pl9 <- pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, percent)) +
  geom_col(fill = "#56B4E9", alpha = 0.8) +
  scale_x_continuous(breaks = 1:9, limits = c(0,10)) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_minimal_hgrid(12)

pl7 | (pl8 / pl9)


