library(tidyverse)
library(magrittr)

raw_data = read_tsv("./_raw/GSE13937_series_matrix.txt",
                    skip = 35,
                    col_names = FALSE)

# Set column and row names and transpose the data
raw_data_flipped = raw_data %>% 
  mutate(X1 = X1 %>% str_replace("!", "")) %>%
  distinct(X1,
           .keep_all = TRUE) %>% 
  column_to_rownames("X1") %>% 
  t() %>% # Transpose, we should probably do this with tidyverse
  as_tibble() %>% 
  column_to_rownames("Sample_title")

raw_data_clean = raw_data_flipped %>% 
  

raw_data_flipped %>% 
  count(Sample_channel_count)