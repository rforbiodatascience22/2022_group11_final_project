# Install required packages -----------------------------------------------


required_packages = c("tidyverse", "ggrepel", "magrittr", "patchwork", "broom",
                      "gridExtra", "survival", "survminer", "cowplot",
                      "DiagrammeRsvg", "rsvg", "DiagrammeR")
for (i in 1:length(required_packages)){
  package = required_packages[i]
  if (!(package %in% as.character(installed.packages()))){
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
source(file = "R/100_kaplan_meier_1.R")
source(file = "R/101_kaplan_meier_2.R")
source(file = "R/104_glm.R")

rmarkdown::render(input = "doc/presentation.Rmd",
                  output_file = "presentation.html")
rmarkdown::render(input = "README.Rmd",
                  output_file = "README.md")
