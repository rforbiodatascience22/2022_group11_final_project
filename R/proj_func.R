# Generate survival probability for Kaplan-Meier plots
generate_survival_probabilities = function(survival_days, shift = FALSE){
  survival_probs = numeric()
  n_observations = length(survival_days)
  for (i in 1:n_observations){
    n_bigger = length(which(survival_days[i] <= survival_days))
    survival_prob = n_bigger / n_observations
    survival_probs = c(survival_probs, survival_prob)
  }
  
  if(shift){
    for(i in 1:length(survival_probs)){
      # Randomize the data a little to avoid duplicates
      survival_probs[i] = survival_probs[i] * runif(n = 1,
                                                    min = 0.999995,
                                                    max = 1.000005)
    }
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

