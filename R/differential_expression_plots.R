# Load libraries ----------------------------------------------------------

library(tidyverse)
library(magrittr)

# Load data ---------------------------------------------------------------

#NOTE: "./data/data_tidy.csv" IS NOW "./data/data_aug_wide.csv".
#YOU SHOULD CHANGE IT IN YOUR READ_CSV FUNCTION
tidy_data = read_csv("./data/data_tidy.csv")

# ADC CT vs NCT  ----------------------------------------------------------

tidy_data %>% 
  select(contains("375"))

tidy_data %>%
  filter(tumor_type == "ADC") %>%
  select(patient:death_due_to_cancer,
         `hsa-mir-21No1`,
         `hsa-mir-223-prec`,
         `hsa-mir-192-2/3No1`,
         `hsa-mir-194-2No1`,
         `hsa-mir-203-precNo1`,
         `hsa-mir-146bNo1`) %>%
  pivot_longer(cols = contains("hsa-mir"),
               names_to = "miRNA",
               values_to = "Expression") %>%
  separate(col = miRNA,
           into = c(NA, NA, "miRNA", NA),
           sep = "-") %>%
  mutate(miRNA = case_when(str_sub(miRNA,
                                   start = -4,
                                   end = -2) == "bNo" ~ str_sub(miRNA,
                                                                end = -5),
                           str_sub(miRNA,
                                   start = -3,
                                   end = -2) == "No" ~ str_sub(miRNA,
                                                                end = -4),
                           TRUE ~ miRNA),
         miRNA = str_c("miR-", miRNA),
         tissue_type = str_to_sentence(tissue_type)) %>%
  ggplot(mapping = aes(x = tissue_type,
                       y = Expression,
                       fill = tissue_type)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.1,
              shape = 21,
              alpha = 0.5,
              show.legend = FALSE) +
  facet_wrap(facets = "miRNA",
             nrow = 2,
             ncol = 3) +
  theme_bw() +
  labs(x = "Tissue type",
       y = "miRNA expression\n(normalized probe intensity)",
       title = "Differentially expressed miRNAs: cancerous (ADC) vs non-cancerous tissue",
       fill = "Tissue type")

ggsave(filename = "./doc/dif_exp_ADC_CT_vs_NCT_2.png",
       plot = last_plot(),
       dpi = 500)

# SCC CT vs NCT -----------------------------------------------------------

tidy_data %>%
  filter(tumor_type == "SCC") %>%
  select(patient:death_due_to_cancer,
         `hsa-mir-21No1`,
         `hsa-mir-375`) %>%
  pivot_longer(cols = contains("hsa-mir"),
               names_to = "miRNA",
               values_to = "Expression") %>%
  separate(col = miRNA,
           into = c(NA, NA, "miRNA", NA),
           sep = "-") %>%
  mutate(miRNA = case_when(str_sub(miRNA,
                                   start = -3,
                                   end = -2) == "No" ~ str_sub(miRNA,
                                                               end = -4),
                           TRUE ~ miRNA),
         miRNA = str_c("miR-", miRNA),
         tissue_type = str_to_sentence(tissue_type)) %>%
  ggplot(mapping = aes(x = tissue_type,
                       y = Expression,
                       fill = tissue_type)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.1,
              shape = 21,
              alpha = 0.5,
              show.legend = FALSE) +
  facet_wrap(facets = "miRNA",
             nrow = 2,
             ncol = 3) +
  theme_bw(base_size = 8) +
  labs(x = "Tissue type",
       y = "miRNA expression\n(normalized probe intensity)",
       title = "Differentially expressed miRNAs: cancerous (SCC) vs non-cancerous tissue",
       fill = "Tissue type")

# ggsave(filename = "./doc/dif_exp_SCC_CT_vs_NCT.png",
#        plot = last_plot(),
#        dpi = 500)

# ADC CT vs SCC CT --------------------------------------------------------

tidy_data %>%
  filter(tissue_type == "cancerous") %>%
  select(patient:death_due_to_cancer,
         `hsa-mir-194-2No1`,
         `hsa-mir-375`) %>%
  pivot_longer(cols = contains("hsa-mir"),
               names_to = "miRNA",
               values_to = "Expression") %>%
  separate(col = miRNA,
           into = c(NA, NA, "miRNA", NA),
           sep = "-") %>%
  mutate(miRNA = case_when(str_sub(miRNA,
                                   start = -3,
                                   end = -2) == "No" ~ str_sub(miRNA,
                                                               end = -4),
                           TRUE ~ miRNA),
         miRNA = str_c("miR-", miRNA),
         tissue_type = str_to_sentence(tissue_type)) %>%
  ggplot(mapping = aes(x = tumor_type,
                       y = Expression,
                       fill = tumor_type)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.1,
              shape = 21,
              alpha = 0.5,
              show.legend = FALSE) +
  facet_wrap(facets = "miRNA",
             nrow = 2,
             ncol = 3) +
  theme_bw(base_size = 8) +
  labs(x = "Tissue type",
       y = "miRNA expression\n(normalized probe intensity)",
       title = "Differentially expressed miRNAs: SCC vs ADC (cancerous tissue)",
       fill = "Tissue type")

# ggsave(filename = "./doc/dif_exp_SCC_vs_ADC.png",
#        plot = last_plot(),
#        dpi = 500)
