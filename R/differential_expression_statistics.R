# Load libraries ----------------------------------------------------------

library("tidyverse")
library("magrittr")
library("broom")
source("./R/proj_func.R")

wide_data = read_csv("./data/data_aug_wide.csv")

# This selects only miRNAs with more than 25% missing samples, as the paper does
valid_miRNAs = wide_data %>%
  select(tumor_type, tissue_type, `hsa-let-7a-1-prec`:ncol(.)) %>% 
  select_if(~sum(is.na(.x)) < nrow(wide_data) * 0.25)

valid_miRNAs_long = valid_miRNAs %>% 
  pivot_longer(cols = starts_with("hsa"),
               names_to = "probe",
               values_to = "expression")

# ADC vs SCC (cancer tissue) ----------------------------------------------

ADC_SCC_dif_exp = diff_exp_stats(filtered_tibble = valid_miRNAs_long %>%
                                        filter(tissue_type == "cancerous"),
                                      group = "tumor_type")

volcano_plot(ttest_tibble = ADC_SCC_dif_exp,
             significance_threshold = 0.05,
             label_logpval_threshold = 4,
             label_logFC_threshold = 0.15,
             plot_title = "SCC vs ADC (cancerous tissues)")

# ADC cancer tissue vs non-cancer tissue ----------------------------------

ADC_CT_NCT_diffexp = diff_exp_stats(filtered_tibble = valid_miRNAs_long %>%
                                           filter(tumor_type == "ADC"),
                                         group = "tissue_type")

volcano_plot(ttest_tibble = ADC_CT_NCT_diffexp,
             significance_threshold = 0.05,
             label_logpval_threshold = 2.5,
             label_logFC_threshold = 0.1)

# SCC cancer tissue vs non-cancer tissue ----------------------------------

SCC_CT_NCT_diffexp = diff_exp_stats(filtered_tibble = valid_miRNAs_long %>%
                                           filter(tumor_type == "SCC"),
                                         group = "tissue_type")

volcano_plot(ttest_tibble = SCC_CT_NCT_diffexp,
             significance_threshold = 0.05,
             label_logpval_threshold = 5,
             label_logFC_threshold = 0.1) +
  xlim(-0.25, 0.25)


# To check if a specific miRNA is differentially expressed between tumor types
# ADC_SCC_significance %>% 
#   filter(BH_pval < 0.05) %>% 
#   select(probe) %>% 
#   filter(str_detect(probe, "155"))
