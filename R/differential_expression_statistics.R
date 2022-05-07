# Load libraries ----------------------------------------------------------

library("tidyverse")
library("magrittr")
library("broom")
source("./R/proj_func.R")

wide_data = read_csv("./data/data_aug_wide.csv")

# This selects only miRNAs with more than 25% missing samples, as the paper does
valid_miRNAs = wide_data %>%
  select(tumor_type, tissue_type, `hsa-let-7a-1-prec`:ncol(.)) %>% 
  select_if(~sum(is.na(.x)) < nrow(wide_data) * 0.25)

valid_miRNAs_long = valid_miRNAs %>% 
  pivot_longer(cols = starts_with("hsa"),
               names_to = "probe",
               values_to = "expression")

valid_miRNAs_nested_cancerous = valid_miRNAs_long %>%
  filter(tissue_type == "cancerous") %>% 
  group_by(probe) %>% 
  nest() %>% 
  ungroup()

valid_miRNAs_tts = valid_miRNAs_nested_cancerous %>% 
  mutate(tt = map(data,
                  ~ t.test(formula = expression ~ tumor_type,
                           data = .x)),
         tt = map(tt,
                  ~ tidy(.x))) %>% 
  unnest(tt) %>%
  rename(ADC_exp = estimate1,
         SCC_exp = estimate2) %>% 
  select(-data)

FDR = 0.1
ADC_SCC_significance = valid_miRNAs_tts %>%
  arrange(p.value) %>% 
  mutate(rank = 1:nrow(.)) %>% 
  mutate(BH_pval = rank / nrow(.) * FDR)  # Adjusted pval by Benjamini-Hochberg method


volcano_plot(ttest_tibble = ADC_SCC_significance,
             significance_threshold = 0.05,
             label_logpval_threshold = 4,
             label_logFC_threshold = 0.15,
             plot_title = "SCC vs ADC (cancerous tissues)")


# To check if a specific miRNA is differentially expressed betwen tumor types
ADC_SCC_significance %>% 
  filter(BH_pval < 0.05) %>% 
  select(probe) %>% 
  filter(str_detect(probe, "155"))
