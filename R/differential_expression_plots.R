# Load libraries ----------------------------------------------------------

library("tidyverse")
library("magrittr")
source("./R/proj_func.R")

# Load data ---------------------------------------------------------------

wide_data = read_csv("./data/data_aug_wide.csv")

# ADC CT vs NCT  ----------------------------------------------------------

diff_exp_boxplot(filtered_tibble = wide_data %>%
                   filter(tumor_type == "ADC"),
                 probes = c("hsa-mir-21No1", "hsa-mir-223-prec",
                            "hsa-mir-192-2/3No1", "hsa-mir-194-2No1", 
                            "hsa-mir-203-precNo1", "hsa-mir-146bNo1"),
                 group = "tissue_type") +
  facet_wrap(facets = "miRNA",
             nrow = 2,
             ncol = 3) +
  labs(x = "Tissue type",
       title = "Differentially expressed miRNAs: cancerous (ADC) vs. 
       non-cancerous tissue",
       fill = "Tissue type")

ggsave(filename = "./result/dif_exp_ADC_CT_vs_NCT_3.png",
       plot = last_plot(),
       dpi = 500)

# SCC CT vs NCT -----------------------------------------------------------

diff_exp_boxplot(filtered_tibble = wide_data %>%
                   filter(tumor_type == "SCC"),
                 probes = c("hsa-mir-21No1", "hsa-mir-375"),
                 group = "tissue_type") +
  facet_wrap(facets = "miRNA",
             nrow = 2,
             ncol = 3) +
  labs(title = "Differentially expressed miRNAs: cancerous (SCC) vs non-cancerous tissue",
       fill = "Tissue type") +
  theme_bw(base_size = 8)

ggsave(filename = "./result/dif_exp_SCC_CT_vs_NCT.png",
       plot = last_plot(),
       dpi = 500)

# ADC CT vs SCC CT --------------------------------------------------------

diff_exp_boxplot(filtered_tibble = wide_data %>%
                   filter(tissue_type == "cancerous"),
                 probes = c("hsa-mir-194-2No1", "hsa-mir-375"),
                 group = "tumor_type") +
  facet_wrap(facets = "miRNA",
             nrow = 2,
             ncol = 3) +
  labs(title = "Differentially expressed miRNAs: SCC vs ADC (cancerous tissue)",
       fill = "Tumor type") +
  theme_bw(base_size = 8)

ggsave(filename = "./result/dif_exp_SCC_vs_ADC.png",
       plot = last_plot(),
       dpi = 500)
