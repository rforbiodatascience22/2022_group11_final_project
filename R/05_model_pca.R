library("tidyverse")
library("patchwork")
library("broom")
library("cowplot")

# Load the files
tidy_data_wide = read_csv("./data/data_aug_wide.csv")
tidy_data_long = read_csv("./data/data_aug_long.csv")

tidy_data_wide <- tidy_data_wide %>%
  filter(tissue_type == "cancerous")

text_size = 18 # for all plots
## PCA
pca_fit <- tidy_data_long %>%
  select(patient, Expression, Mature_MiRNA, tissue_type) %>%
  filter(tissue_type == "cancerous") %>%
  group_by(Mature_MiRNA, patient) %>%
  mutate(Mean = mean(Expression)) %>%
  select(-Expression) %>%
  distinct() %>%
  pivot_wider(names_from = Mature_MiRNA, values_from = Mean) %>%
  select_if(~ !any(is.na(.))) %>%
  ungroup() %>%
  select(-patient, -tissue_type) %>%
  scale() %>%
  prcomp(scale = TRUE)

pl7 <- pca_fit %>%
  augment(tidy_data_wide) %>%
  ggplot(aes(.fittedPC1,
             .fittedPC2,
             color = cohort)) +
  geom_point(size = 3,
             mapping = aes(shape = cohort)) +
  labs(x = "PC1",
       y = "PC2",
       title = "PCA analysis") +
  theme_half_open(12) +
  background_grid()

## Rotation matrix
# Define arrow style for plotting
arrow_style <- arrow(
  angle = 20, ends = "first", type = "closed", length = grid::unit(8, "pt")
)

# Plot the matrix
pl8 <- pca_fit %>%
  tidy(matrix = "rotation") %>%
  pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value") %>%
  ggplot(aes(PC1, PC2)) +
  geom_segment(xend = 0,
               yend = 0,
               arrow = arrow_style) +
  geom_text(
    aes(label = column),
    hjust = 1,
    nudge_x = -0.02,
    color = "#904C2F",
    size = text_size / 3,
    check_overlap = T
  ) +
  labs(title = "Rotation matrix"
  ) +
  xlim(-0.3, .3) + ylim(-.3, 0.3) +
  theme_minimal_grid(12)

## Plot the principal components
pl9 <- pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, percent)) +
  geom_col(fill = "#56B4E9", alpha = 0.8) +
  labs(title = "Contribution of components"
  ) +
  scale_x_continuous(breaks = 1:9, limits = c(0,10)) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_minimal_hgrid(12)

# Combine the plots into one figure
pca_plots <- (pl7 | (pl8 / pl9)) & theme(text = element_text(size = text_size),
                                         axis.text = element_text(size = text_size),
                                         legend.text = element_text(size = text_size))
pca_plots +
  plot_annotation(theme = theme(plot.title = element_text(size = text_size * 2.15,
                                                          hjust = 0.5,
                                                          margin = margin(0,
                                                                          0,
                                                                          10,
                                                                          0))))

ggsave(filename = "./results/PCA1.png",
       plot = last_plot(),
       width = 4000,
       height = 2000,
       units = "px")
