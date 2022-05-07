volcano_plot = function(ttest_tibble, significance_threshold,
                        label_logpval_threshold, label_logFC_threshold,
                        plot_title = ""){
  library("ggrepel")
  v_plot = ttest_tibble %>% 
    mutate(FC = ADC_exp / SCC_exp,
           logFC = log2(FC),
           logP = -log10(p.value),
           significance = case_when(BH_pval < significance_threshold ~
                                      "Significant",
                                    TRUE ~ "Not significant"),
           significance = significance %>% 
             factor(labels = c("Significant", "Not significant")),
           conditional_label = case_when(logP > label_logpval_threshold &
                                           abs(logFC) > label_logFC_threshold ~ 
                                           probe,
                                         TRUE ~ "")) %>%
    ggplot(mapping = aes(x = logFC,
                         y = logP,
                         fill = significance,
                         label = conditional_label)) +
    geom_point(shape = 21,
               alpha = 0.7) +
    # geom_text(hjust=0,
    #           vjust=-1) +
    geom_text_repel(force = 150) +
    geom_vline(xintercept = label_logFC_threshold,
               linetype = "dashed",
               alpha = 0.2) +
    geom_vline(xintercept = -label_logFC_threshold,
               linetype = "dashed",
               alpha = 0.2) +
    geom_hline(yintercept = label_logpval_threshold,
               linetype = "dashed",
               alpha = 0.2) +
    scale_fill_manual(values = c("orange", "purple")) +
    theme_bw() +
    labs(title = plot_title,
         x = "log(FC)",
         y = "-log(p value)",
         fill = "Significance")
  return(v_plot)
}

diff_exp_boxplot = function(filtered_tibble, probes, group, plot_title = ""){
  plot = filtered_tibble %>%  # Must be wide version
    select(patient:death_due_to_cancer,
           probes) %>%
    pivot_longer(cols = contains("hsa-mir"),
                 names_to = "miRNA",
                 values_to = "Expression") %>%
    separate(col = miRNA,
             into = c(NA, NA, "miRNA", NA),
             sep = "-") %>%
    mutate(miRNA = case_when(str_sub(miRNA,
                                     start = -3,
                                     end = -2) == "No" ~ str_sub(miRNA,
                                                                 end = -4),
                             TRUE ~ miRNA),
           miRNA = str_c("miR-", miRNA),
           tissue_type = str_to_sentence(tissue_type)) %>%
    ggplot(mapping = aes(x = case_when(group == "tumor_type" ~ tumor_type,
                                       group == "tissue_type" ~ tissue_type),
                         y = Expression,
                         fill = case_when(group == "tumor_type" ~ tumor_type,
                                          group == "tissue_type" ~ tissue_type))) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(outlier.alpha = 0) +
    geom_jitter(width = 0.1,
                shape = 21,
                alpha = 0.5,
                show.legend = FALSE) +
    theme_bw() +
    labs(x = group,
         y = "miRNA expression\n(normalized probe intensity)",
         title = plot_title,
         fill = group)
  return(plot)
}
