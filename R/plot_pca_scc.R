library(tidyverse)
library(patchwork)
library(broom)
library(cowplot)
library(ggpubr)

# Load the files
tidy_data_wide = read_csv("./data/data_aug_wide.csv")
tidy_data_long = read_csv("./data/data_aug_long.csv")

# Filter for squamous cell carcinoma only
tidy_data_long <- tidy_data_long %>%
  filter(tumor_type == "SCC")
toDelete <- seq(0, length(tidy_data_wide), 2)
tidy_data_wide <- tidy_data_wide[-toDelete,] %>%
  filter(tumor_type == "SCC")

## PCA
pca_fit <- tidy_data_long %>%
  select(patient, Expression, Mature_MiRNA) %>%
  group_by(Mature_MiRNA, patient) %>%
  mutate(Mean = mean(Expression)) %>%
  select(-Expression) %>%
  distinct() %>%
  pivot_wider(names_from = Mature_MiRNA, values_from = Mean) %>%
  select_if(~ !any(is.na(.))) %>%
  ungroup() %>%
  select(-patient) %>%
  scale() %>%
  prcomp(scale = TRUE)

pl7 <- pca_fit %>%
  augment(tidy_data_wide) %>%
  ggplot(aes(.fittedPC1, .fittedPC2, color = nodal_involvement)) +
  geom_point(size = 3) +
  labs(x = "PC1",
       y = "PC2",
       title = "PCA analysis") +
  theme_half_open(12) +
  background_grid()

pl7

# ## Rotation matrix
# # Define arrow style for plotting
# arrow_style <- arrow(
#   angle = 20, ends = "first", type = "closed", length = grid::unit(8, "pt")
# )
# 
# # Plot the matrix
# pl8 <- pca_fit %>%
#   tidy(matrix = "rotation") %>%
#   pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value") %>%
#   ggplot(aes(PC1, PC2)) +
#   geom_segment(xend = 0, yend = 0, arrow = arrow_style) +
#   geom_text(
#     aes(label = column),
#     hjust = 1, nudge_x = -0.02, 
#     color = "#904C2F", size = 2,
#     check_overlap = T
#   ) +
#   labs(title = "Rotation matrix"
#   ) +
#   xlim(-0.3, .3) + ylim(-.3, 0.3) +
#   theme_minimal_grid(12)
# 
# ## Plot the principal components 
# pl9 <- pca_fit %>%
#   tidy(matrix = "eigenvalues") %>%
#   ggplot(aes(PC, percent)) +
#   geom_col(fill = "#56B4E9", alpha = 0.8) +
#   labs(title = "Contribution of components"
#   ) +
#   scale_x_continuous(breaks = 1:9, limits = c(0,10)) +
#   scale_y_continuous(
#     labels = scales::percent_format(),
#     expand = expansion(mult = c(0, 0.01))
#   ) +
#   theme_minimal_hgrid(12)
# 
# pl7 | (pl8 / pl9)

# ggsave(filename = "./results/PCA_SCC.png",
#        plot = last_plot(),
#        dpi = 1000,
#        width = 15,
#        height = 10)
