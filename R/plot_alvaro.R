library(tidyverse)
library(patchwork)
library(magrittr)


data_tidy_long<- read_csv("data/data_aug_long.csv") 

data_tidy_wide = read_csv("./data/data_aug_wide.csv") 

view(data_tidy_filtered)

data_tidy

data_tidy_filtered <- data_tidy_filtered %>% 
  drop_na()



#Attempt to do Figure 2 
#A) A) Altered expression between cancerous tissue (NCT) and non-cancerous tissue (NCT) in ADC patients 
#is shown. Mir-21, mir-223, mir-192, and mir-194are over-expressed in tumors while mir-203 is under-expressed in tumors. 


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


#Next: create graph!! Heatmap??
data <- data %>% 
  mutate(ratio = non_cancerous / cancerous) %>% 
  filter(ratio > 1) %>% 
  arrange(desc(ratio)) %>% 
  relocate(probe, .before = non_cancerous)


















### ------------------------------ THIS DOES NOT WORK VERY WELL --------
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
              values_from = expression) %>% 
  group_by(probe) %>% View()

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



data_tidy <- data_tidy_filtered %>% 
  select(-"Mature_MiRNA") %>% 
  pivot_wider(names_from = "Probe_name",
              values_from = "Expression")




#Removing the "Probe name" column so that we only have the "Mature MiRNA"
#data_tidy_filtered <- data_tidy_filtered %>% 
#  select(-"Probe_name")



#Warning: the following lines are not well finished.


#Pivoting back to the different MiRNAs into columns
# Apparently, the mature MiRNAs are duplicated or triplicated in the
#microarray, so as we do not have only 1 entry per observation, you can't really repivot it. At least I have
#not been able to do it and make it work. I suggest working with data_tidy_filtered after we decide what to do
#with the duplicates and triplicates, but I would understand if you think there are too many observations
#(I don't think they are that many now)
# data_tidy_repivoted <- data_tidy_filtered %>% 
#   pivot_wider(names_from = "Mature_MiRNA",
#               values_from = "Expression" )
# 

#This is to check the IDs/Number of mature MiRNAs that we should have 
# mature_data = read_tsv("./_raw/A-GEOD-13415.adf.txt",
#                        skip = 15,
#                        col_names = TRUE)

