# Insert linear model here, Jake
library(tidyverse)
library(patchwork)
library(broom)
library(cowplot)

# Load the files
tidy_data_wide = read_csv("./data/data_aug_wide.csv")
tidy_data_long = read_csv("./data/data_aug_long.csv")

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