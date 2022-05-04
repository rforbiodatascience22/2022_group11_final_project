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


### Filtering only the relevant MiRNAs of the paper ###
## Loading the data of the microarray probes ##
data_tidy <- read_csv("./data/data_tidy.csv")

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



