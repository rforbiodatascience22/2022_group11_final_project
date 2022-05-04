library(tidyverse)
library(patchwork)
library(magrittr)
source("./R/proj_func.R")

data_tidy_filtered <- read_csv("data/data_tidy_filtered.csv") %>% 
  as_tibble() 

data_tidy <- read_csv("data/data_tidy.csv") %>% 
  as_tibble() 

view(data_tidy_filtered)

data_tidy

data_tidy_filtered <- data_tidy_filtered %>% 
  drop_na()



#Attempt to do Figure 2 
#A) A) Altered expression between cancerous tissue (NCT) and non-cancerous tissue (NCT) in ADC patients 
#is shown. Mir-21, mir-223, mir-192, and mir-194are over-expressed in tumors while mir-203 is under-expressed in tumors. 


First I have to select the samples with adenocarcinoma.
The x axis is mir-21, mir-223, mir-203, mir-192, mir-194 
The y axis is Expression of cancerous tissue vs non-cancerous tissue (in log2)

Stratify by Country (UMD/Canada)

data_tidy %>% 
  select("tumor_type",
         "cohort",
         "hsa-mir-21No1",
         "hsa-mir-223-prec",
         "hsa-mir-203-precNo1",
         "hsa-mir-192-2/3No1") %>% 
  pivot_longer(cols = contains("hsa"),
               names_to = "probe",
               values_to = "expression") %>% 
  ggplot(mapping = aes(x = probe,
                       y = expression)) +
  geom_col(mapping = aes(fill = cohort), 
           position="dodge")
  

non_grouped <- data_tidy %>%
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
  pivot_longer(cols = contains("hsa"),
               names_to = "probe",
               values_to = "expression") %>%
  pivot_wider(names_from = tissue_type,
              values_from = expression)

%>% 
  group_by(probe) 

%>% 
  mutate(differential_expresion = "cancerous" - "non-cancerous") %>% View()


#Me he quedado aqui arriba justo. Tengo que restar la expresion de canceroso a la de no canceroso. 

  #filter(tumor_type == "SCC") %>% 
  ggplot(mapping = aes(x = probe,
                       y = expression)) +
  geom_col(mapping = aes(fill = tissue_type), 
           position="dodge") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_wrap(vars(cohort),
      ncol = 2,
      scales = "free_y")





ggplot(data = data_tidy_filtered,
       mapping = aes(x = smoking,
                     y = survival_days)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(mapping = aes(fill = nodal_involvement),
             shape = 21,
             color = "black") 
#facet_wrap(vars(species),
#    ncol = 2,
#   scales = "free_y")





#What I will be doing:
#  - Download txt file from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GPL8835&id=30413&db=GeoDb_blob35 to get the miRNA names
#- I should probably end up with something like this https://www.ebi.ac.uk/arrayexpress/files/A-GEOD-13415/A-GEOD-13415.adf.txt


