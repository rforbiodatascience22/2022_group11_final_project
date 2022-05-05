# Load libraries ----------------------------------------------------------
library("tidyverse")
library("magrittr")

# Objective: Add new variables to your data

# Define functions --------------------------------------------------------
#source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
#my_data_clean <- read_tsv(file = "data/02_my_data_clean.tsv")

data_tidy <- read_csv("./data/data_tidy.csv")

# Wrangle data ------------------------------------------------------------
#my_data_clean_aug <- my_data_clean # %>% ...

# Grouping based on staging (in situ, localized, regional, and distant)
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



### Filtering only the relevant MiRNAs of the paper ###
## Loading the data of the microarray probes ##

raw_probes <- read_tsv("./_raw/A-GEOD-8835.adf.txt",
                       skip = 14,
                       col_names = TRUE)

raw_comments <- read_tsv("./_raw/A-GEOD-8835_comments.txt",
                         col_names = TRUE)

raw_probes <- bind_cols(raw_probes,
                        raw_comments)

#This is the control probes - we do NOT want these -
raw_controls <- raw_probes %>% 
  filter(`Comment[SPOT_ID]` == 'Control') 

#This is the raw_probes without all the controls. From now on, we no longer need the 4 variables created above.
#We will only work with probes_data
probes_data <- anti_join(raw_probes,
                         raw_controls)

#Now I am only selecting the human data and the data from mature MiRNAs, and 
#then Removing unnecessary columns so that we only have data for mature MiRNAs
probes_data <- probes_data %>% 
  filter(`Reporter Group [organism]` == "Homo sapiens" &
           `Comment[Contains_Mature_MiRNA]` == "yes"  &
           str_detect(`Reporter Name`, 'hsa')) %>% 
  select(-c(`Comment[SPOT_ID]`,`Reporter Group [organism]`,`Comment[Contains_Mature_MiRNA]`))

#Renaming the two columns
probes_data <- probes_data %>% 
  rename("Probe_name" = "Reporter Name")  %>% 
  rename("Mature_MiRNA" = "Reporter Database Entry [mirbase]")

#Pivoting the data_tidy miRNAs to be able to filter them later on 
data_tidy_pivoted <- data_tidy %>% 
  pivot_longer(cols = "ath-MIR156aNo1" : "series_matrix_table_end", 
               names_to = "Probe_name", 
               values_to = "Expression")

#Filtering the rows from data_tidy that have a match in probes_data  (I need to keep both columns in probes_data)
#from probes_data so that I change the names of the MiRNA next
data_tidy_filtered <- right_join(data_tidy_pivoted,
                                 probes_data) %>% 
  drop_na() %>% 
  relocate(Expression, .after = last_col())

write_csv(data_tidy_filtered, 
          file = "./data/data_tidy_filtered.csv")


data_tidy_filtered <- read_csv("data/data_tidy_filtered.csv") %>% 
  as_tibble() 


#Pivoting back to cartesian data tidy
data_tidy <- data_tidy_filtered %>% 
  select(-"Mature_MiRNA") %>% 
  pivot_wider(names_from = "Probe_name",
              values_from = "Expression")


# Write data --------------------------------------------------------------
# write_tsv(x = my_data_clean_aug,
#           file = "data/03_my_data_clean_aug.tsv")

write_csv(data_tidy, 
          file = "./data/data_tidy.csv")