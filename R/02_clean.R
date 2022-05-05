# Load libraries ----------------------------------------------------------
library("tidyverse")
library("magrittr")

# Objective: Remove invalid data, e.g. if you have amino acid sequence data, 
# remove non-valid sequences containing X or other non-standard amino acid characters 
# or fix columns, e.g. dates or when two labels are the same, but spelled differently

# Define functions --------------------------------------------------------
#source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
raw_data_flipped <- read_csv("./data/data_transposed.csv")
# my_data <- read_tsv(file = "data/01_my_data.tsv")

# Wrangle data ------------------------------------------------------------

# Remove unnecessary columns, most of these have all the same value or are
# completely irrelevant
data_concise = raw_data_flipped %>% 
  select(-c(Sample_status:Sample_channel_count,
            Sample_organism_ch1,
            hybridization_date,
            histology,
            Sample_molecule_ch1:Sample_platform_id,
            Sample_contact_name:Sample_contact_country,
            Sample_supplementary_file:ID_REF,
            Blank))


data_tidy = data_concise %>% 
  separate(col = Sample_title,
           into = c(NA, "patient")) %>% 
  mutate(patient = str_sub(string = patient,
                           end = -2)) %>% 
  separate(col = Sample_source_name_ch1,
           sep = " ",
           into = c("tumor_type", NA, "tissue_type")) %>% 
  mutate(cohort = str_sub(string = cohort,
                          start = 9)) %>%
  mutate(barretts_esophagus = str_sub(string = barretts_esophagus,
                                      start = 22)) %>%
  mutate(chemoradiation_therapy = str_sub(string = chemoradiation_therapy,
                                          start = 25)) %>% 
  mutate(alcohol_consumption = str_sub(string = alcohol_consumption,
                                       start = 21),
         alcohol_consumption = na_if(alcohol_consumption,
                                     y = "")) %>%
  mutate(smoking = str_sub(string = smoking,
                           start = 10),
         smoking = na_if(smoking,
                         y = "")) %>%
  mutate(ptnm_stage = str_sub(string = ptnm_stage,
                              start = 13)) %>%
  mutate(nodal_involvement = str_sub(string = nodal_involvement,
                                     start = 20)) %>%
  mutate(survival_days = str_sub(string = survival_days,
                                 start = 16),
         survival_days = as.integer(survival_days)) %>% 
  mutate(death_due_to_cancer = str_sub(string = death_due_to_cancer,
                                       start = 22))





# Write data --------------------------------------------------------------
# write_tsv(x = my_data_clean,
#           file = "data/02_my_data_clean.tsv")

write_csv(data_tidy, 
          file = "./data/data_tidy.csv")
