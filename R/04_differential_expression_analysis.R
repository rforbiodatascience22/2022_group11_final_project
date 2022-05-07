# Load libraries ----------------------------------------------------------

library("tidyverse")
library("magrittr")
library("broom")
source("./R/proj_func.R")

# Prepare data ------------------------------------------------------------

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

logpval_threshold = 4
logFC_threshold = 0.15

volcano_plot(ttest_tibble = ADC_SCC_dif_exp,
             significance_threshold = 0.05,
             label_logpval_threshold = logpval_threshold,
             label_logFC_threshold = logFC_threshold,
             plot_title = "SCC vs ADC (cancerous tissues)")

ggsave(filename = "./results/ADC_SCC_volc_plot.png",
       height = 13,
       width = 13,
       units = "cm",
       dpi = 500)

plot_probes = get_plot_candidates(ttest_tibble = ADC_SCC_dif_exp,
                                  logpval_threshold = logpval_threshold,
                                  logFC_threshold = logFC_threshold)

diff_exp_boxplot(filtered_tibble = wide_data %>%
                   filter(tissue_type == "cancerous"),
                 group = "tumor_type",
                 probes = plot_probes,
                 facet_rows = 1,
                 facet_scale = "free_y") +
  labs(title = "SCC vs ADC (cancerous tissues)")
  theme_bw(base_size = 13.5) +
  theme(legend.position = "none")

ggsave(filename = "./results/ADC_SCC_diffexp_plot.png",
       height = 10,
       width = 25,
       units = "cm",
       dpi = 500)

paper_probes = ADC_SCC_dif_exp %>%
  filter(str_detect(string = probe,
                    pattern = "194|375")) %>% 
  pull(probe)

diff_exp_boxplot(filtered_tibble = wide_data %>%
                   filter(tissue_type == "cancerous"),
                 group = "tumor_type",
                 probes = paper_probes,
                 facet_rows = 1,
                 facet_scale = "free_y") +
  labs(title = "SCC vs ADC (cancerous tissues)",
       subtitle = "Differentially expressed miRNAs from Mathé et al.") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

ggsave(filename = "./results/ADC_SCC_diffexp_plot_paper.png",
       height = 10,
       width = 14,
       units = "cm",
       dpi = 500)

# ADC cancer tissue vs non-cancer tissue ----------------------------------

ADC_CT_NCT_diffexp = diff_exp_stats(filtered_tibble = valid_miRNAs_long %>%
                                           filter(tumor_type == "ADC"),
                                         group = "tissue_type")

logpval_threshold = 2.5
logFC_threshold = 0.1

volcano_plot(ttest_tibble = ADC_CT_NCT_diffexp,
             significance_threshold = 0.05,
             label_logpval_threshold = logpval_threshold,
             label_logFC_threshold = logFC_threshold,
             plot_title = "ADC cancerous vs. non cancerous tissue") +
  xlim(-0.4, 0.2)  +
  ylim(0, 3.6)

ggsave(filename = "./results/ADC_CT_NCT_volc_plot.png",
       height = 13,
       width = 13,
       units = "cm",
       dpi = 500)

plot_probes = get_plot_candidates(ttest_tibble = ADC_CT_NCT_diffexp,
                                  logpval_threshold = logpval_threshold,
                                  logFC_threshold = logFC_threshold)

diff_exp_boxplot(filtered_tibble = wide_data %>%
                   filter(tumor_type == "ADC"),
                 group = "tissue_type",
                 probes = plot_probes,
                 facet_rows = 1) +
  labs(title = "ADC cancerous vs. non cancerous tissue") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

ggsave(filename = "./results/ADC_CT_NCT_diffexp_plot.png",
       height = 10,
       width = 25,
       units = "cm",
       dpi = 500)

paper_probes = ADC_CT_NCT_diffexp %>%
  filter(str_detect(string = probe,
                    pattern = "-21[N-]|223|192|194|203")) %>% 
  pull(probe)

diff_exp_boxplot(filtered_tibble = wide_data %>%
                   filter(tumor_type == "ADC"),
                 group = "tissue_type",
                 probes = paper_probes,
                 facet_rows = 1) +
  labs(title = "ADC cancerous vs. non cancerous tissue",
       subtitle = "Differentially expressed miRNAs from Mathé et al.") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

ggsave(filename = "./results/ADC_CT_NCT_diffexp_plot_paper.png",
       height = 10,
       width = 25,
       units = "cm",
       dpi = 500)

# SCC cancer tissue vs non-cancer tissue ----------------------------------

SCC_CT_NCT_diffexp = diff_exp_stats(filtered_tibble = valid_miRNAs_long %>%
                                           filter(tumor_type == "SCC"),
                                         group = "tissue_type")

logpval_threshold = 5
logFC_threshold = 0.1

volcano_plot(ttest_tibble = SCC_CT_NCT_diffexp,
             significance_threshold = 0.05,
             label_logpval_threshold = logpval_threshold,
             label_logFC_threshold = logFC_threshold,
             plot_title = "ADC cancerous vs. non cancerous tissue") +
  xlim(-0.25, 0.25)

ggsave(filename = "./results/SCC_CT_NCT_volc_plot.png",
       height = 13,
       width = 13,
       units = "cm",
       dpi = 500)

plot_probes = get_plot_candidates(ttest_tibble = SCC_CT_NCT_diffexp,
                                  logpval_threshold = logpval_threshold,
                                  logFC_threshold = logFC_threshold)

diff_exp_boxplot(filtered_tibble = wide_data %>%
                   filter(tumor_type == "SCC"),
                 group = "tissue_type",
                 probes = plot_probes,
                 facet_rows = 2) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  labs(title = "SCC cancerous vs. non cancerous tissue")

ggsave(filename = "./results/SCC_CT_NCT_diffexp_plot.png",
       height = 18,
       width = 25,
       units = "cm",
       dpi = 500)

paper_probes = SCC_CT_NCT_diffexp %>%
  filter(str_detect(string = probe,
                    pattern = "-21[N-]|375")) %>% 
  pull(probe)

diff_exp_boxplot(filtered_tibble = wide_data %>%
                   filter(tumor_type == "SCC"),
                 group = "tissue_type",
                 probes = paper_probes,
                 facet_rows = 1) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none") +
  labs(title = "SCC cancerous vs. non cancerous tissue",
       subtitle = "Differentially expressed miRNAs from Mathé et al.")

ggsave(filename = "./results/SCC_CT_NCT_diffexp_plot_paper.png",
       height = 10,
       width = 14,
       units = "cm",
       dpi = 500)

# Specific manual checks --------------------------------------------------

# To check if a specific miRNA is differentially expressed between tumor types
# ADC_SCC_dif_exp %>%
#   filter(BH_pval < 0.05) %>%
#   select(probe) %>%
#   filter(str_detect(probe, "155"))

# To check if a specific miRNA is differentially expressed between tissue types
# ADC_CT_NCT_diffexp %>% 
#   filter(BH_pval < 0.05) %>%
#   select(probe) %>%
#   filter(str_detect(probe, "155"))

# SCC_CT_NCT_diffexp %>% 
#   filter(BH_pval < 0.05) %>%
#   select(probe) %>%
#   filter(str_detect(probe, "155"))
