library(tidyverse)
library(patchwork)
library(broom)
library(cowplot)

tidy_data_wide = read_csv("./data/data_aug_wide.csv")
tidy_data_long = read_csv("./data/data_aug_long.csv")

# toDelete <- seq(0, length(tidy_data_wide), 2)
# data <- tidy_data_wide[-toDelete,]
#   
# 
# pl1 <- data %>%
#   ggplot(aes(x = survival_days, y = smoking, fill = death_due_to_cancer)) +
#   geom_boxplot(alpha = 0.5) +
#   theme_classic(base_family = "Avenir",
#                 base_size = 12) +
#   theme(legend.position = "none") + 
#   labs(x = "Survival days",
#        y = "Smoking",
#        title = str_c("Survival days and cause of death of patients under ",
#                      "different coditions"
#                     )
#        )
# 
# pl2 <- data %>%
#   ggplot(aes(x = survival_days, y = alcohol_consumption, fill = death_due_to_cancer)) +
#   geom_boxplot(alpha = 0.5) +
#   theme_classic(base_family = "Avenir",
#                 base_size = 12) +
#   theme(legend.position = "none") + 
#   labs(x = "Survival days",
#        y = "Alcohol consumption"
#        )
# 
# pl3 <- data %>%
#   ggplot(aes(x = survival_days, y = chemoradiation_therapy, fill = death_due_to_cancer)) +
#   geom_boxplot(alpha = 0.5) +
#   theme_classic(base_family = "Avenir",
#                 base_size = 12) +
#   theme(legend.position = "bottom") + 
#   labs(x = "Survivay days",
#        y = "Chemo therapy",
#        fill = "Death due to cancer",
#        )
# 
# pl1 / pl2 / pl3
# 
# pl4 <- data %>%
#   ggplot(aes(x = survival_days,
#              y = smoking,
#              fill = death_due_to_cancer)) +
#   geom_boxplot(alpha = 0.5) +
#   labs(x = "Survival days",
#        y = "Smoking") +
#   theme_classic(base_family = "Avenir",
#                 base_size = 12) +
#   theme(legend.position = "none")
# 
# pl5 <- data %>%
#   ggplot(aes(x = survival_days,
#              y = alcohol_consumption,
#              fill = death_due_to_cancer)) +
#   geom_boxplot(alpha = 0.5) +
#   labs(x = "Survival days",
#        y = "Alcohol consumption") +
#   theme_classic(base_family = "Avenir",
#                 base_size = 12) +
#   theme(legend.position = "none")
# 
# pl6 <- data %>%
#   group_by(death_due_to_cancer) %>%
#   summarise(`Chemo therapy` = sum(chemoradiation_therapy == "Yes"), 
#             Smoking = sum(smoking == "Yes"),
#             `Alcohol consumption` = sum(alcohol_consumption == "Yes"),
#             `Barretts esophagus` = sum(barretts_esophagus == "Yes"),
#             `Nodal involvement` = sum(nodal_involvement == "Yes")) %>%
#   pivot_longer(!death_due_to_cancer) %>%
#   ggplot(aes(x = value, y = name, fill = death_due_to_cancer)) +
#   geom_col(position = "dodge", alpha = 0.5) +
#   labs(x = "Patient count",
#        y = "Conditions",
#        fill = "Death due to cancer") +
#   theme_classic(base_family = "Avenir",
#               base_size = 12) + 
#   theme(legend.position = "bottom")
# 
# pl6 / (pl4 + pl5)

pca_plot <- function(crit){
  # PCA
  pca_fit <- tidy_data_long %>% 
    select(patient, Expression, Mature_MiRNA) %>%
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
  
  data_aug <- tidy_data_wide[-toDelete,]
  
  pl7 <- pca_fit %>%
    augment(data_aug) %>%
    ggplot(aes_string(".fittedPC1", ".fittedPC2", color = crit)) + 
    geom_point(size = 1.5) +
    labs(x = "PC1",
         y = "PC2",
         title = "PCA analysis") +
    theme_half_open(12) + 
    background_grid()
  
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
    labs(title = "Rotation matrix"
    ) +
    xlim(-0.3, .3) + ylim(-.3, 0.3) +
    theme_minimal_grid(12)
  
  pl9 <- pca_fit %>%
    tidy(matrix = "eigenvalues") %>%
    ggplot(aes(PC, percent)) +
    geom_col(fill = "#56B4E9", alpha = 0.8) +
    labs(title = "Contribution of components"
    ) +
    scale_x_continuous(breaks = 1:9, limits = c(0,10)) +
    scale_y_continuous(
      labels = scales::percent_format(),
      expand = expansion(mult = c(0, 0.01))
    ) +
    theme_minimal_hgrid(12)
  
  pl7 | (pl8 / pl9)
}

pca_plot(crit = "tumor_type")

pca_plot1 <- function(crit){
  pca_fit <- tidy_data_wide %>% 
    select_if(~ !any(is.na(.))) %>%
    select(matches("hsa")) %>%
    scale() %>%
    prcomp(scale = TRUE)
  
  pl7 <- pca_fit %>%
    augment(tidy_data_wide) %>%
    ggplot(aes_string(".fittedPC1", ".fittedPC2", color = crit)) + 
    geom_point(size = 1.5) +
    labs(x = "PC1",
         y = "PC2",
         title = "PCA analysis") +
    theme_half_open(12) + 
    background_grid()
  pl7
}

pca_plot1(crit = "SEER_stage")

ggsave(filename = "./results/PCA.png",
       plot = last_plot(),
       dpi = 1000,
       width = 15,
       height = 10)
# Linear-model

data1 <- tidy_data_long %>% 
  select(tissue_type, Expression, Mature_MiRNA) %>%
  group_by(Mature_MiRNA, tissue_type) %>%
  mutate(Mean = mean(Expression)) %>%
  select(-Expression) %>%
  distinct() %>%
  pivot_wider(names_from = tissue_type, values_from = Mean) %>%
  ungroup() %>%
  drop_na()

mod <- glm(data1$cancerous ~ data1$`non-cancerous`)

data1 <- mod %>% 
  augment(data1) %>%
  mutate(`Differently expressed` = case_when(.resid < -2 ~ "Low",
                        -2 < .resid & .resid < 2 ~ "Neutral",
                        2 <= .resid ~ "High"))

pl10 <- data1 %>%
  ggplot(aes(`non-cancerous`, cancerous, label = Mature_MiRNA, color=`Differently expressed`)) +
  geom_point(size = 1.5) + 
  geom_smooth(method='glm', color='black', linetype = 'longdash') +
  geom_text(aes(label=ifelse(.resid > 2 | .resid < -2, as.character(Mature_MiRNA), '')), 
            hjust=0, vjust=0) +
  labs(x = "Non-cancerous",
       y = "Cancerous",
       title = "Mean miRNA expression of non-cancerous vs cancerous and SCC vs ADC",
       subtitle = "Diffentially expressed miRNAs are labeled") +
  theme_half_open(12) + 
  theme(legend.position = "none")

data2 <- tidy_data_long %>% 
  select(tumor_type, Expression, Mature_MiRNA) %>%
  group_by(Mature_MiRNA, tumor_type) %>%
  mutate(Mean = mean(Expression)) %>%
  select(-Expression) %>%
  distinct() %>%
  pivot_wider(names_from = tumor_type, values_from = Mean) %>%
  ungroup() %>%
  drop_na()

mod <- glm(data2$ADC ~ data2$SCC)

data2 <- mod %>% 
  augment(data2) %>%
  mutate(`Differently expressed` = case_when(.resid < -2 ~ "Low",
                                             -2 < .resid & .resid < 2 ~ "Neutral",
                                             2 <= .resid ~ "High"))

pl11 <- data2 %>%
  ggplot(aes(SCC, ADC, label = Mature_MiRNA, color=`Differently expressed`)) +
  geom_point(size = 1.5) + 
  geom_smooth(method='glm', color='black', linetype = 'longdash') +
  geom_text(aes(label=ifelse(.resid > 2 | .resid < -2, as.character(Mature_MiRNA), '')), 
            hjust=0, vjust=0) +
  theme_half_open(12)

pl10 + pl11

ggsave(filename = "./results/lm1.png",
       plot = last_plot(),
       dpi = 1000,
       width = 15,
       height = 10)

data3 <- tidy_data_long %>% 
  select(Expression, Mature_MiRNA, tissue_type, tumor_type) %>%
  filter(tumor_type == "ADC") %>%
  group_by(Mature_MiRNA, tissue_type) %>%
  mutate(Mean = mean(Expression)) %>%
  select(-Expression) %>%
  distinct() %>%
  pivot_wider(names_from = tissue_type, values_from = Mean) %>%
  ungroup() %>%
  drop_na()

mod <- glm(data3$cancerous ~ data3$`non-cancerous`)

data3 <- mod %>% 
  augment(data3) %>%
  mutate(`Differently expressed` = case_when(.resid < -2 ~ "Low",
                                             -2 < .resid & .resid < 2 ~ "Neutral",
                                             2 <= .resid ~ "High"))


pl12 <- data3 %>%
  ggplot(aes(`non-cancerous`, cancerous, label = Mature_MiRNA, color=`Differently expressed`)) +
  geom_point(size = 1.5) + 
  geom_smooth(method='glm', color='black', linetype = 'longdash') +
  geom_text(aes(label=ifelse(.resid > 2 | .resid < -2, as.character(Mature_MiRNA), '')), 
            hjust=0, vjust=0) +
  labs(x = "Non-cancerous",
       y = "Cancerous",
       title = "ADC") +
  theme_half_open(12) +
  theme(legend.position = "none")

data4 <- tidy_data_long %>% 
  select(Expression, Mature_MiRNA, tissue_type, tumor_type) %>%
  filter(tumor_type == "SCC") %>%
  group_by(Mature_MiRNA, tissue_type) %>%
  mutate(Mean = mean(Expression)) %>%
  select(-Expression) %>%
  distinct() %>%
  pivot_wider(names_from = tissue_type, values_from = Mean) %>%
  ungroup() %>%
  drop_na()

mod <- glm(data4$cancerous ~ data4$`non-cancerous`)

data4 <- mod %>% 
  augment(data4) %>%
  mutate(`Differently expressed` = case_when(.resid < -2 ~ "Low",
                                             -2 < .resid & .resid < 2 ~ "Neutral",
                                             2 <= .resid ~ "High"))
pl13 <- data4 %>%
  ggplot(aes(`non-cancerous`, cancerous, label = Mature_MiRNA, color=`Differently expressed`)) +
  geom_point(size = 1.5) + 
  geom_smooth(method='glm', color='black', linetype = 'longdash') +
  geom_text(aes(label=ifelse(.resid > 2 | .resid < -2, as.character(Mature_MiRNA), '')), 
            hjust=0, vjust=0) +
  labs(x = "Non-cancerous",
       y = "Cancerous",
       title = "SCC") +
  theme_half_open(12)

pl12 + pl13

ggsave(filename = "./results/lm2.png",
       plot = last_plot(),
       dpi = 1000,
       width = 15,
       height = 10)

data5 <- tidy_data_long %>% 
  select(Expression, Mature_MiRNA, tissue_type, tumor_type) %>%
  filter(tissue_type == "cancerous") %>%
  group_by(Mature_MiRNA, tumor_type) %>%
  mutate(Mean = mean(Expression)) %>%
  select(-Expression) %>%
  distinct() %>%
  pivot_wider(names_from = tumor_type, values_from = Mean) %>%
  ungroup() %>%
  drop_na()

mod <- glm(data5$SCC ~ data5$ADC)

data5 <- mod %>% 
  augment(data5) %>%
  mutate(`Differently expressed` = case_when(.resid < -2 ~ "Low",
                                             -2 < .resid & .resid < 2 ~ "Neutral",
                                             2 <= .resid ~ "High"))


pl14 <- data5 %>%
  ggplot(aes(ADC, SCC, label = Mature_MiRNA, color=`Differently expressed`)) +
  geom_point(size = 1.5) + 
  geom_smooth(method='glm', color='black', linetype = 'longdash') +
  geom_text(aes(label=ifelse(.resid > 2 | .resid < -2, as.character(Mature_MiRNA), '')), 
            hjust=0, vjust=0) +
  labs(title = "Cancerous") +
  theme_half_open(12) +
  theme(legend.position = "none")

data6 <- tidy_data_long %>% 
  select(Expression, Mature_MiRNA, tissue_type, tumor_type) %>%
  filter(tissue_type == "non-cancerous") %>%
  group_by(Mature_MiRNA, tumor_type) %>%
  mutate(Mean = mean(Expression)) %>%
  select(-Expression) %>%
  distinct() %>%
  pivot_wider(names_from = tumor_type, values_from = Mean) %>%
  ungroup() %>%
  drop_na()

mod <- glm(data6$SCC ~ data6$ADC)

data6 <- mod %>% 
  augment(data6) %>%
  mutate(`Differently expressed` = case_when(.resid < -2 ~ "Low",
                                             -2 < .resid & .resid < 2 ~ "Neutral",
                                             2 <= .resid ~ "High"))
pl15 <- data6 %>%
  ggplot(aes(ADC, SCC, label = Mature_MiRNA, color=`Differently expressed`)) +
  geom_point(size = 1.5) + 
  geom_smooth(method='glm', color='black', linetype = 'longdash') +
  geom_text(aes(label=ifelse(.resid > 2 | .resid < -2, as.character(Mature_MiRNA), '')), 
            hjust=0, vjust=0) +
  labs(title = "Non-cancerous") +
  theme_half_open(12)

pl14 + pl15

ggsave(filename = "./results/lm3.png",
       plot = last_plot(),
       dpi = 1000,
       width = 15,
       height = 10)
