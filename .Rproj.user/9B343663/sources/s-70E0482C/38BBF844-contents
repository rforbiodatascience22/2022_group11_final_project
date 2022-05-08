# Load libraries ----------------------------------------------------------

library("tidyverse")
library("magrittr")

# Load raw data -----------------------------------------------------------

raw_data = read_tsv("./_raw/GSE13937_series_matrix.txt",
                    skip = 35,
                    col_names = FALSE)

# Set column and row names and transpose the data
raw_data_flipped = raw_data %>% 
  mutate(X1 = X1 %>% str_replace("!", "")) %>%
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


##Loading the data of the microarray probes #

raw_probes <- read_tsv("./_raw/A-GEOD-8835.adf.txt",
                       skip = 14,
                       col_names = TRUE)

raw_comments <- read_tsv("./_raw/A-GEOD-8835_comments.txt",
                         col_names = TRUE)

probes_data <- bind_cols(raw_probes,
                        raw_comments)


# Write data --------------------------------------------------------------
write_csv(raw_data_flipped, 
          file = "./data/data_load.csv")

write_csv(probes_data,
          file = "./data/probes_data_load.csv")
