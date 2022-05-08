# Generates volcano plot from differential expression analysis
volcano_plot = function(ttest_tibble, significance_threshold,
                        label_logpval_threshold, label_logFC_threshold,
                        plot_title = ""){
  library("ggrepel")
  v_plot = ttest_tibble %>% 
    mutate(significance = case_when(BH_pval < significance_threshold ~
                                      "Significant",
                                    TRUE ~ "Not significant"),
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
    geom_text_repel(force = 110) +
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
    theme(legend.position = "bottom",
          text = element_text(family = "Helvetica")) +
    labs(title = plot_title,
         x = "log(FC)",
         y = "-log(p value)",
         fill = "Significance")
  return(v_plot)
}

# Generates fancy box plot from expression data
diff_exp_boxplot = function(filtered_tibble, probes, group, facet_rows = NULL,
                            facet_cols = NULL, facet_scale = "free_y"){
  probes_clean = read_csv("data/probes_data_clean.csv") %>% 
    rename(probe = Probe_name,
           miRNA = Mature_MiRNA)
  
  plot = filtered_tibble %>%  # Must be wide version
    select(patient:death_due_to_cancer,
           probes) %>%
    pivot_longer(cols = contains("hsa-mir"),
                 names_to = "probe",
                 values_to = "Expression") %>%
    inner_join(y = probes_clean,
               by = "probe") %>% 
    mutate(miRNA = str_sub(miRNA,
                           start = 5),
           tissue_type = str_to_sentence(tissue_type),
           tissue_type = case_when(tissue_type == "Non-cancerous" ~ "Healthy",
                                   TRUE ~ tissue_type)) %>%
    ggplot(mapping = aes(x = case_when(group == "tumor_type" ~ tumor_type,
                                       group == "tissue_type" ~ tissue_type),
                         y = Expression,
                         fill = case_when(group == "tumor_type" ~ tumor_type,
                                          group == "tissue_type" ~ 
                                            tissue_type))) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(outlier.alpha = 0) +
    geom_jitter(width = 0.1,
                shape = 21,
                alpha = 0.5,
                show.legend = FALSE) +
    theme_bw() +
    labs(x = case_when(group == "tumor_type" ~ "Tumor type",
                       group == "tissue_type" ~ "Tissue type"),
         y = "miRNA expression\n(normalized probe intensity)",
         fill = case_when(group == "tumor_type" ~ "Tumor type",
                          group == "tissue_type" ~ "Tissue type")) +
    facet_wrap(facets = "miRNA",
               nrow = facet_rows,
               ncol = facet_cols,
               scales = facet_scale) +
    theme(legend.position = "none",
          text = element_text(family = "Helvetica"))
  return(plot)
}

# Performs the differential expression analysis as described in the paper
diff_exp_stats = function(filtered_tibble, group, FDR = 0.1){
  nested = filtered_tibble %>% 
    group_by(probe) %>% 
    nest() %>% 
    ungroup()
  
  t_tests = nested %>% 
    mutate(tt = map(data,
                    ~ t.test(formula = expression ~ 
                               case_when(group == "tumor_type" ~ tumor_type,
                                         group == "tissue_type" ~ tissue_type),
                             data = .x)),
           tt = map(tt, ~ tidy(.x))) %>% 
    unnest(tt) %>%
    rename(ADC_exp = estimate1,
           SCC_exp = estimate2) %>% 
    select(-data)
  
  results = t_tests %>%
    arrange(p.value) %>% 
    mutate(rank = 1:nrow(.),
           BH_pval = rank / nrow(.) * FDR, # Adj pval Benjamini-Hochberg method
           FC = ADC_exp / SCC_exp,
           logFC = log2(FC),
           logP = -log10(p.value))  
  
  return(results)
}

# Returns a vector with the probe names of differentially expressed miRNAS,
# based on given criteria
get_plot_candidates = function(ttest_tibble, logpval_threshold,
                               logFC_threshold){
  probes = ttest_tibble %>% 
    filter(abs(logFC) > logFC_threshold & logP > logpval_threshold) %>% 
    pull(probe)
  return(probes)
}

# Generates survival probability for Kaplan-Meier plots
# OLD: Not tidyverse and dropped in favor of using the survival package. Still
# here so kaplan_meier.R still works (dead end, though).
generate_survival_probabilities = function(survival_days, shift = FALSE){
  survival_probs = numeric()
  n_observations = length(survival_days)
  for (i in 1:n_observations){
    n_bigger = length(which(survival_days[i] <= survival_days))
    survival_prob = n_bigger / n_observations
    survival_probs = c(survival_probs, survival_prob)
  }
  
  if(shift){
    # Randomize the data a little to avoid duplicates
    survival_probs = survival_probs * runif(n = length(survival_probs),
                                            min = 0.999995,
                                            max = 1.000005)

    survival_shifted = numeric(length(survival_probs))
    
    sorted_probs = sort(survival_probs,
                        decreasing = TRUE)
    
    for(i in 1:length(survival_probs)){
      if(i == 1){
        prev_survival = 1
      }
      
      ith = which(sorted_probs[i] == survival_probs)
      
      survival_shifted[ith] = prev_survival
      
      prev_survival = sorted_probs[i]
    }
    
    survival_probs = survival_shifted
  }
  
  return(survival_probs)
}
