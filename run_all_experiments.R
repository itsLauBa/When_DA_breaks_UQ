# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

status_coverage_ex <- system2("Rscript", file.path("coverage_example_section4.R"))
if (status_coverage_ex != 0) {
  stop("Coverage example failed")
}

status_classical <- system2("Rscript", file.path("classical_experiment.R"))
if (status_classical != 0) {
  stop("Classical experiment failed")
}

status_high_dim <- system2("Rscript", file.path("high_dim_experiment.R"))
if (status_high_dim != 0) {
  stop("High-dimensional experiment failed")
}

status_plots <- system2("Rscript", file.path("plot_results.R"))
if (status_plots != 0) {
  stop("Plot generation failed")
}
