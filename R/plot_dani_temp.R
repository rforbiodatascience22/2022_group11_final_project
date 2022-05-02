library(tidyverse)
library(magrittr)
source("./R/proj_func.R")

tidy_data = read_csv("./data/data_tidy.csv")

colnames(tidy_data)

data = tidy_data %>% 
  select(patient:death_due_to_cancer,
         `hsa-mir-21No1`,
         `hsa-mir-223-prec`,
         `hsa-mir-192-2/3No1`,
         `hsa-mir-194-2No1`,
         `hsa-mir-203-precNo1`,
         `hsa-mir-375`
         )  


data %>%
  filter(tumor_type == "SCC") %>% 
  ggplot(mapping = aes(x = tissue_type,
                       y = `hsa-mir-21No1`,
                       fill = tissue_type)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.1,
              shape = 21,
              show.legend = FALSE) +
  theme_bw()

# ggsave(filename = "./doc/dif-exp.png",
#        plot = last_plot(),
#        dpi = 500)


# Between cancers
data %>%
  filter(tissue_type == "cancerous") %>% 
  ggplot(mapping = aes(x = `hsa-mir-375`,
                       fill = tumor_type)) +
  geom_density(alpha = 0.5) +
  xlim(5, 18)

## Cancerous vs. non-cancerous
data %>%
  filter(tumor_type == "SCC") %>% 
  ggplot(mapping = aes(x = `hsa-mir-21No1`,
                       fill = tissue_type)) +
  geom_density(alpha = 0.5) +
  xlim(5, 15)
