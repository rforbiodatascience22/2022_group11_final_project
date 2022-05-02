library(tidyverse)
library(magrittr)
source("./R/proj_func.R")

data_tidy_filtered <- read_csv("data/data_tidy_filtered.csv") %>% 
  as_tibble() %>% 
  drop_na()

view(data_tidy_filtered)


data_tidy_filtered <- data_tidy_filtered %>% 
  drop_na()



#Attempt to do Figure 2 
#A) A) Altered expression between cancerous tissue (NCT) and non-cancerous tissue (NCT) in ADC patients 
#is shown. Mir-21, mir-223, mir-192, and mir-194are over-expressed in tumors while mir-203 is under-expressed in tumors. 







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


