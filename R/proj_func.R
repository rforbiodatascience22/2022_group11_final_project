# Generate survival probability for Kaplan-Meier plots
generate_survival_probabilities = function(survival_days){
  # ??? Consult Leon on how to do this in tidyverse
  survival_probs = numeric()
  n_observations = length(survival_days)
  for (i in 1:n_observations){
    n_bigger = length(which(survival_days[i] <= survival_days))
    survival_prob = n_bigger / n_observations
    survival_probs = c(survival_probs, survival_prob)
  }
  return(survival_probs)
}

