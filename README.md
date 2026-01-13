# Internal Validation Simulation
Version 2 of the internal validation simulation project.

`simulation.R` defines a data-generating process for a binary outcome using
specified coefficients, sets the intercept to achieve 15% prevalence, fits a
logistic model on a simulated sample of 631 observations, and evaluates
performance (AUC, calibration slope, Brier score, and MAPE) on a large
simulated dataset.
