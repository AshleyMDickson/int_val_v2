# Internal Validation Simulation: Version 2

This project conducts an extensive simulation study to evaluate and compare internal validation methods for clinical prediction models (logistic regression). It focuses on how accurately different techniques estimate "true" external performance across multiple metrics.

## Project Overview

The simulation explores the performance of three common internal validation methods—**Optimism-Corrected Bootstrapping**, **10-fold Cross-Validation**, and **Sample Splitting**—against the "Ground Truth" calculated from a massive independent dataset.

### Performance Metrics Evaluated:
- **Calibration Slope (CS)**: Measures the agreement between predicted and observed risks (Target = 1.0).
- **C-statistic / AUC**: Measures the discriminative ability of the model.
- **Brier Score**: Measures the overall accuracy of probabilistic predictions.
- **Mean Absolute Prediction Error (MAPE)**: Measures the average magnitude of prediction error.

## Simulation Design

- **Data Generating Process (DGP)**: 10 predictors with specified coefficients and a tuned intercept to achieve ~15% outcome prevalence. The true AUC is approximately 0.75.
- **Sample Sizes**: The simulation evaluates three scenarios based on the target shrinkage of 0.9:
  - **Small ($N=518$)**: 50% of the recommended size.
  - **Recommended ($N=1036$)**: Calculated to achieve a 0.9 calibration slope.
  - **Large ($N=1554$)**: 150% of the recommended size.
- **Iterations**: Each scenario is repeated 200 times.
- **Validation**: "True" performance is measured on a separate dataset of 500,000 observations.

## Key Findings (Summary)

Based on the `summary_statistics.csv` results:

1.  **Bootstrap Superiority**: The Bootstrap method consistently provides the most unbiased estimates of the true external performance. For the recommended sample size ($N=1036$), the Bootstrap mean Calibration Slope (~0.913) almost perfectly matches the External True Slope (~0.913).
2.  **The Reliability Paradox**: While 10-fold Cross-Validation is relatively unbiased on average, it often shows a "reliability paradox"—its point-by-point correlation with the truth in individual datasets can be noisier or less stable than the Bootstrap.
3.  **Sample Splitting Inefficiency**: Sample splitting (50/50) consistently underestimates performance and shows significantly higher variance (Standard Deviation) compared to the other two methods. This is expected as it only utilizes half the data for model development.
4.  **Convergence**: As sample size increases, the gap between internal estimates and the external truth narrows for all methods, and the standard deviation of estimates decreases significantly.

## Visualizations

The script generates 8 primary visualizations:
- **`scatter_{metric}.png`**: Displays the correlation between internal estimates and external truth. A 1:1 dashed line indicates perfect agreement.
- **`box_{metric}.png`**: Shows the distribution of estimates across the 200 simulations, with jittered points representing individual results.

## Files
- `simulation.R`: Main script for the DGP, validation loops, and plotting.
- `simulation_results.csv`: Raw data from all simulations.
- `summary_statistics.csv`: Mean and SD for all methods and metrics.
- `sim_helpers.R`: Utility functions for convergence checks and plotting.

---
*Created for the Internal Validation V2 project.*
