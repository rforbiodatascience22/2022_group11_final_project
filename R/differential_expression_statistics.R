# Load libraries ----------------------------------------------------------

library("tidyverse")
library("magrittr")
# library("rstatix")

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
  unnest(tt)



valid_miRNAs %>% 
  filter(tissue_type == "cancerous") %>% 
  select(-tissue_type) %>%
  t_test(formula = `hsa-let-7a-1-prec` ~ tumor_type)


