# Load libraries ----------------------------------------------------------
library("tidyverse")
library("magrittr")

# Load data ---------------------------------------------------------------
#my_data_clean <- read_tsv(file = "data/02_my_data_clean.tsv")

data_tidy <- read_csv("./data/data_clean.csv")
probes_data <- read_csv("./data/probes_data_clean.csv")

# Wrangle data ------------------------------------------------------------
#my_data_clean_aug <- my_data_clean # %>% ...

# Grouping data_tidy based on staging (in situ, localized, regional, and distant)
data_tidy <- data_tidy %>% 
  mutate(SEER_stage = 
           case_when(ptnm_stage == "0" ~ "In situ",
                     nodal_involvement == "No" & ptnm_stage != "0" ~ "Localized",
                     nodal_involvement == "Yes" & ptnm_stage != "IV" ~ "Regional",
                     ptnm_stage == "IV" ~ "Distant")) %>%
  mutate(SEER_stage = factor(SEER_stage,
                             levels = c("In situ",
                                        "Localized",
                                        "Regional",
                                        "Distant"))) %>%
  relocate(SEER_stage, .after = ptnm_stage)




#Pivoting the data_tidy miRNAs to be able to filter them later on 
data_tidy_pivoted <- data_tidy %>% 
  pivot_longer(cols = "ath-MIR156aNo1" : "series_matrix_table_end", 
               names_to = "Probe_name", 
               values_to = "Expression")

#Filtering the rows from data_tidy that have a match in probes_data  
#(I need to keep both columns in probes_data)
#from probes_data so that I change the names of the MiRNA next
data_tidy_filtered <- right_join(data_tidy_pivoted,
                                 probes_data) %>% 
  drop_na() %>% 
  relocate(Expression, .after = last_col())


#Pivoting back to cartesian data tidy
data_tidy <- data_tidy_filtered %>% 
  select(-"Mature_MiRNA") %>% 
  pivot_wider(names_from = "Probe_name",
              values_from = "Expression")


# Write data --------------------------------------------------------------
# write_tsv(x = my_data_clean_aug,
#           file = "data/03_my_data_clean_aug.tsv")
write_csv(data_tidy_filtered, 
          file = "./data/data_aug_long.csv")

write_csv(data_tidy, 
          file = "./data/data_aug_wide.csv")
