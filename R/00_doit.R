# Install required packages -----------------------------------------------

required_packages = c("tidyverse", "ggrepel")
for (i in 1:length(needed_packages)){
  package = needed_packages[i]
  if (!(package %in% installed.packages() %>% as.character())){
    install.packages(package)
  }
}

# Run all scripts ---------------------------------------------------------
source(file = "R/01_load.R")
source(file = "R/02_clean.R")
source(file = "R/03_augment.R")
source(file = "R/04_model_kaplan.R")
source(file = "R/05_model_pca.R")
source(file = "R/06_model_diffexp.R")
source(file = "R/07_model_diffexp_plot.R")
source(file = "R/08_flowchart.R")
rmarkdown::render(input = "doc/presentation.Rmd",
                  output_file = "presentation.html")
rmarkdown::render(input = "README.Rmd",
                  output_file = "README.md")
