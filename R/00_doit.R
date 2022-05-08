# Run all scripts ---------------------------------------------------------
source(file = "R/01_load.R")
source(file = "R/02_clean.R")
source(file = "R/03_augment.R")
source(file = "R/04_differential_expression_analysis.R")
source(file = "R/99_diffexp_plot_composit.R")
rmarkdown::render(input = "doc/presentation.Rmd",
                  output_file = "presentation.html")
