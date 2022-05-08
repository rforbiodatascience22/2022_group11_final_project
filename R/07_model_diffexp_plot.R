library("patchwork")
source("./R/06_model_diffexp.R")

vplots = ADC_SCC_vplot +
  labs(x = NULL,
       title = "ADC vs. SCC\n(cancerous tissues)") +
  ADC_CT_NCT_vplot +
  labs(y = NULL,
       title = "Adenocarcinoma\n(cancerous vs. non-cancerous)") +
  SCC_CT_NCT_vplot +
  labs(x = NULL,
       y = NULL,
       title = "Squamous cell carcinoma\n(cancerous vs. non-cancerous)")

vplots +
  plot_layout(guides = "collect") &
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(filename = "./results/presentation_vplots.png",
       dpi = 400,
       width = 28,
       height = 13,
       units = "cm")
