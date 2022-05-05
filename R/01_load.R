# Load libraries ----------------------------------------------------------

library("tidyverse")
library("magrittr")

# Objective: Collapse data to a single file or convert .xlsx to .tsv, here we could 
# imagine having an .xlsx-file with multiple sheets, from which we create a single .tsv

# Load raw data -----------------------------------------------------------

raw_data = read_tsv("./_raw/GSE13937_series_matrix.txt",
                    skip = 35,
                    col_names = FALSE)

# Set column and row names and transpose the data
raw_data_flipped = raw_data %>% 
  mutate(X1 = X1 %>% str_replace("!", "")) %>% 
  # If there's time, try to find better way of doing this 
  mutate(X1 = case_when(
    str_sub(X2, end = 6) == "cohort" ~ "cohort",
    str_sub(X2, end = 18) == "hybridization date" ~ "hybridization_date",
    str_sub(X2, end = 11) == "tissue type" ~ "tissue_type",
    str_sub(X2, end = 9) == "histology" ~ "histology",
    str_sub(X2, end = 19) == "barrett's esophagus" ~ "barretts_esophagus",
    str_sub(X2, end = 22) == "chemoradiation therapy" ~ 
      "chemoradiation_therapy",
    str_sub(X2, end = 19) == "alcohol consumption" ~ "alcohol_consumption",
    str_sub(X2, end = 7) == "smoking" ~ "smoking",
    str_sub(X2, end = 10) == "ptnm stage" ~ "ptnm_stage",
    str_sub(X2, end = 17) == "nodal involvement" ~ "nodal_involvement",
    str_sub(X2, end = 13) == "survival days" ~ "survival_days",
    str_sub(X2, end = 19) == "death due to cancer" ~ "death_due_to_cancer",
    TRUE ~ X1)) %>%
  column_to_rownames("X1") %>% 
  t() %>%
  as_tibble()



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

#This is the raw_probes without all the controls. 
#From now on, we no longer need the 4 variables created above.
#We will only work with probes_data
probes_data <- anti_join(raw_probes,
                         raw_controls)


# Write data --------------------------------------------------------------
write_csv(raw_data_flipped, 
          file = "./data/data_transposed.csv")

write_csv(probes_data,
          file = "./data/probes_data.csv")
