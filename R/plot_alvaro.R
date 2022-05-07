library(tidyverse)
library(patchwork)
library(magrittr)
library(patchwork)

data_tidy_long<- read_csv("data/data_aug_long.csv") 

data_tidy_wide = read_csv("./data/data_aug_wide.csv") 

#This was an attempted plot that did not come to fruition due to lack of the data used in the
#paper 


#Selecting the most important probes and getting the expression column for cancerous
#and non-cancerous
data_selected <- data_tidy_wide %>%
  select(1:14,
         "hsa-mir-21No1",
         "hsa-mir-021-prec-17No1",
         "hsa-mir-223-prec",
         "hsa-mir-146bNo1",
         "hsa-mir-146-prec",
         "hsa-mir-181a-precNo1",
         "hsa-mir-181b-precNo1",
         "hsa-mir-103-prec-5=103-1",
         "hsa-mir-107-prec-10",
         "hsa-let-7c-prec",
         "hsa-mir-203-precNo1",
         "hsa-mir-205-prec") %>%
  select(-Sample_geo_accession) %>%  #It gives issues with the mapping 
  pivot_longer(cols = contains("hsa"),
               names_to = "probe",
               values_to = "expression") %>%
  pivot_wider(names_from = tissue_type,
              values_from = expression) %>%
  rename("non_cancerous" = "non-cancerous")


#Use of mapping function to remove all the NAs and compress the dataframe
data <-  data_selected  %>% 
  nest(-patient, -probe) %>% 
  mutate(data = map(data, ~ map_dfc(., na.omit))) %>% 
  unnest() %>% 
  distinct() 

#Getting the ratio of expression 
data <- data %>% 
  mutate(ratio = non_cancerous / cancerous) %>% 
  filter(ratio > 1) %>% 
  #arrange(desc(ratio)) %>% 
  relocate(probe, .before = non_cancerous)


#Next: create graph!! Heatmap??
# ggplot(data = data,
#        mapping = aes(x= probe, 
#                      y= patient)) +
#   geom_tile(mapping = aes(fill = ratio)) +

  

