library("patchwork")
source("./R/04_differential_expression_analysis.R")

vplots = ADC_SCC_vplot +
  labs(x = NULL) +
  ADC_CT_NCT_vplot +
  labs(y = NULL) +
  SCC_CT_NCT_vplot +
  labs(x = NULL,
       y = NULL)

vplots +
  plot_layout(guides = "collect") &
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

ggsave(filename = "./results/presentation_vplots.png",
       dpi = 400,
       width = 35,
       height = 15,
       units = "cm")
