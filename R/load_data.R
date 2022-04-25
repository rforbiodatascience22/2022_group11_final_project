library(tidyverse)
library(magrittr)

raw_data = read_tsv("./_raw/GSE13937_series_matrix.txt",
                    skip = 35,
                    col_names = FALSE)

# Set column and row names and transpose the data
raw_data_flipped = raw_data %>% 
  mutate(X1 = X1 %>% str_replace("!", "")) %>% 
  # ??? Better way of doing this? I've tried mutate with case_when (unsuccessfully) 
  mutate(X1 = case_when(str_sub(X2, end = 6) == "cohort" ~ "cohort",
                        str_sub(X2, end = 18) == "hybridization date" ~"hybridization_date", # ??? respect max 80 chars
                        str_sub(X2, end = 11) == "tissue type" ~ "tissue_type",
                        str_sub(X2, end = 9) == "histology" ~ "histology",
                        str_sub(X2, end = 19) == "barrett's esophagus" ~ "barretts_esophagus",
                        str_sub(X2, end = 22) == "chemoradiation therapy" ~ "chemoradiation_therapy",
                        str_sub(X2, end = 19) == "alcohol consumption" ~ "alcohol_consumption",
                        str_sub(X2, end = 7) == "smoking" ~ "smoking",
                        str_sub(X2, end = 10) == "ptnm stage" ~ "ptnm_stage",
                        str_sub(X2, end = 17) == "nodal involvement" ~ "nodal_involvement",
                        str_sub(X2, end = 13) == "survival days" ~ "survival_days",
                        str_sub(X2, end = 19) == "death due to cancer" ~ "death_due_to_cancer",
                        TRUE ~ X1)) %>%
  column_to_rownames("X1") %>% 
  t() %>% # ??? Transpose, is it ok to do this in tidyverse?
  as_tibble()
  
write_csv(raw_data_flipped, 
          file = "./data/data_transposed.csv")


# Remove unnecessary columns, most of these have all the same value or are
# completely irrelevant
data_concise = raw_data_flipped %>% 
  select(-c(Sample_status:Sample_channel_count,
            Sample_organism_ch1,
            hybridization_date,
            histology,
            Sample_molecule_ch1:Sample_platform_id,
            Sample_contact_name:Sample_contact_country,
            Sample_supplementary_file:ID_REF))


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

#Removing the probes coming from Mus musculus
data_tidy <- data_tidy %>% 
  select(-contains("mmu"))

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

write_csv(data_tidy, 
          file = "./data/data_tidy.csv")

