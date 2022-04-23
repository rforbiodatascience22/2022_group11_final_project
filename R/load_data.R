library(tidyverse)
library(magrittr)

raw_data = read_tsv("./_raw/GSE13937_series_matrix.txt",
                    skip = 35,
                    col_names = FALSE)

# Set column and row names and transpose the data
raw_data_flipped = raw_data %>% 
  mutate(X1 = X1 %>% str_replace("!", ""))
  # ??? Better way of doing this? I've tried mutate with case_when (unsuccesfully) 
  mutate(X1 = case_when(str_sub(X2, end = 6) == "cohort" ~ "cohort",
                        str_sub(X2, end = 18) == "hybridization date" ~"hybridization date", # ??? respect max 80 chars
                        str_sub(X2, end = 11) == "tissue type" ~ "tissue type",
                        str_sub(X2, end = 9) == "histology" ~ "histology",
                        str_sub(X2, end = 19) == "barrett's esophagus" ~ "barrett's esophagus",
                        str_sub(X2, end = 22) == "chemoradiation therapy" ~ "chemoradiation therapy",
                        str_sub(X2, end = 19) == "alcohol consumption" ~ "alcohol consumption",
                        str_sub(X2, end = 7) == "smoking" ~ "smoking",
                        str_sub(X2, end = 10) == "ptnm stage" ~ "ptnm stage",
                        str_sub(X2, end = 17) == "nodal involvement" ~ "nodal involvement",
                        str_sub(X2, end = 13) == "survival days" ~ "survival days",
                        str_sub(X2, end = 19) == "death due to cancer" ~ "death due to cancer",
                        TRUE ~ X1)) %>%
  column_to_rownames("X1") %>% 
  t() %>% # ??? Transpose, there must be a tidyverse way of doing this
  as_tibble() %>%
  column_to_rownames("Sample_title")
  
write_csv(raw_data_flipped, 
          file = "./data/data_clean.csv")


# TODO: move this part of the wrangling into a new  
  
# Remove unnecessary columns, most of these have all the same value
raw_data_concise = raw_data_flipped %>% 
  select(-c(Sample_status,
            Sample_submission_date,
            Sample_last_update_date,
            Sample_type,  # Always RNA
            Sample_channel_count,
            Sample_organism_ch1,
            Sample_molecule_ch1,
            Sample_extract_protocol_ch1,
            Sample_label_ch1,
            Sample_label_protocol_ch1,
            Sample_taxid_ch1,
            Sample_hyb_protocol,
            Sample_scan_protocol,
            Sample_description,
            Sample_data_processing,
            Sample_platform_id,
            Sample_contact_name:Sample_contact_country)) %>% View()  # TODO: Add the rest of the irrelevant columns


