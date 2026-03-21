# Import libraries
library(pROC)
#library(rms)
library(samplesizedev)
library(foreach)
library(doParallel)
library(ggplot2)
library(tidyverse)

# # Set-up parameters and such
beta <- c(0.40, 0.35, -0.35, 0.30, -0.25, 0.20, 0.15, 0.10, -0.6, 0.05) # tuned to get AUC ~ 0.75 with 10 predictors, and prevalence ~ 15% (with alpha = -3)
k <- length(beta)
phi <- 0.15
#set.seed(1729)

get_alpha <- function(target_prev, beta, k) {
  N <- 100000
  predictors <- matrix(rnorm(N*k), nrow = N, ncol = k)
  lp_fixed <- predictors %*% beta
  root <- function(a) mean(plogis(a + lp_fixed)) - target_prev
  alpha <- uniroot(root, c(-10, 10))$root
  return(alpha)
}
alpha <- get_alpha(phi, beta, k)
print(paste0("Alpha: ", alpha)) # <- -2.041136

dgp <- function(n, k, return_auc = FALSE) {
  X <- matrix(rnorm(n*k), nrow = n, ncol = k)
  
  linear_predictor <- alpha + (X %*% beta)
  probability <- as.numeric(1/(1+exp(-1*linear_predictor)))
  outcome <- rbinom(n, size = 1, prob = probability)
  
  df <- data.frame(outcome, X)
  
  if (return_auc == TRUE) {
    auc_true <- auc(outcome, probability, quiet = TRUE)
    return(list(df = df, auc_true = auc_true))
  } else {
    return(df)
  }
}

test_data <- dgp(500000, k, return_auc = TRUE)
auc_val <- test_data$auc_true
print(paste0("True AUC: ", auc_val))
print(paste0("Marginal prevalence: ", mean(test_data$df$outcome)))
#samp_size_calc <- suppressWarnings(samplesizedev(outcome = "Binary", S = 0.9, phi = phi, c = auc_val, p = k))
#sample_size <- samp_size_calc$sim
sample_size <- 1036 # calculated from samplesizedev, with 10 predictors, phi = 0.15, AUC = 0.75, and shrinkage = 0.9. This is the sample size needed to get a shrinkage of 0.9 in the final model, which is what we want to use for this simulation. We can also try with a different sample sizes.
print(paste0("Required sample size = ", sample_size))

# Helper function to calculate all performance metrics
calc_metrics <- function(outcome, lp) {
  prob <- plogis(lp)
  # Calibration slope
  cs <- unname(coef(glm(outcome ~ lp, family = "binomial"))[2])
  # AUC
  auc_val <- as.numeric(auc(outcome, prob, quiet = TRUE))
  # Brier Score
  brier <- mean((outcome - prob)^2)
  # MAPE
  mape <- mean(abs(outcome - prob))
  
  return(c(cs = cs, auc = auc_val, brier = brier, mape = mape))
}

# True model performance
# test_data must be supplied as a dataframe; note that dgp() returns a list with df and auc_true, so we need to use test_data$df when calling this function
model_truth <- function(df, test_data) {
  model <- glm(outcome ~ ., family = "binomial", data = df)
  ext_lp <- predict(model, newdata = test_data, type = "link") # test_data drawn once outside this function
  calc_metrics(test_data$outcome, ext_lp)
}

boot_corr_metrics <- function(df, B = 200) {
  n <- nrow(df)
  # Apparent performance
  model_app <- glm(outcome ~ ., family = "binomial", data = df)
  lp_app <- predict(model_app, newdata = df, type = "link")
  perf_app <- calc_metrics(df$outcome, lp_app)
  
  optimism <- matrix(NA, nrow = B, ncol = length(perf_app))
  colnames(optimism) <- names(perf_app)
  
  for(i in 1:B) {
    idx <- sample(n, replace = TRUE)
    df_boot <- df[idx, ]
    model_boot <- glm(outcome ~ ., family = "binomial", data = df_boot)
    
    # Performance on bootstrap sample
    lp_boot_boot <- predict(model_boot, newdata = df_boot, type = "link")
    perf_boot_boot <- calc_metrics(df_boot$outcome, lp_boot_boot)
    
    # Performance on original sample
    lp_boot_orig <- predict(model_boot, newdata = df, type = "link")
    perf_boot_orig <- calc_metrics(df$outcome, lp_boot_orig)
    
    optimism[i, ] <- perf_boot_boot - perf_boot_orig
  }
  
  perf_corr <- perf_app - colMeans(optimism)
  return(perf_corr)
}

sample_split_metrics <- function(df, prop_train = 0.5) {
  n <- nrow(df)
  idx <- sample.int(n)
  n_train <- floor(prop_train * n)
  train_idx <- idx[1:n_train]
  test_idx  <- idx[(n_train + 1):n]
  
  train_df <- df[train_idx, , drop = FALSE]
  test_df  <- df[test_idx,  , drop = FALSE]
  
  model <- glm(outcome ~ ., family = "binomial", data = train_df)
  lp_test <- predict(model, newdata = test_df, type = "link")
  calc_metrics(test_df$outcome, lp_test)
}

cv_metrics <- function(df, K = 10) {
  n <- nrow(df)
  folds <- sample(rep(1:K, length.out = n))
  lp_oof <- rep(NA_real_, n) # pooled out of fold linear preds
  
  for (kfold in 1:K) {
    train_df <- df[folds != kfold, , drop = FALSE]
    test_df  <- df[folds == kfold, , drop = FALSE]
    model <- glm(outcome ~ ., family = "binomial", data = train_df)
    lp_oof[folds == kfold] <- predict(model, newdata = test_df, type = "link")
  }
  calc_metrics(df$outcome, lp_oof)
}

single_comparison <- function(sample_size, i) {
  set.seed(i)
  df <- dgp(sample_size, k)
  test_df <- dgp(100000, k)
  
  perf_true  <- model_truth(df, test_df)
  perf_boot  <- boot_corr_metrics(df, B = 200)
  perf_split <- sample_split_metrics(df, prop_train = 0.5)
  perf_cv    <- cv_metrics(df, K = 10)

  # Combine and name
  res <- c(
    setNames(perf_true, paste0(names(perf_true), "_true")),
    setNames(perf_boot, paste0(names(perf_boot), "_boot")),
    setNames(perf_split, paste0(names(perf_split), "_split")),
    setNames(perf_cv, paste0(names(perf_cv), "_cv"))
  )
  return(res)
}

simulation <- function(sample_size, nsim) { # n = number of simualtion repetitions
  cores <- detectCores() - 1
  cluster <- makeCluster(cores)
  registerDoParallel(cluster)
  print(paste0("Running simulation on ", cores, " cores..."))
  
  simulation_data <- foreach(
    i = 1:nsim,
    .combine = rbind,
    .packages = c("rms", "pROC"),
    .export = c("single_comparison", "dgp", "model_truth", "cv_metrics", 
                "sample_split_metrics", "boot_corr_metrics", "calc_metrics",
                "alpha", "beta", "k")
  ) %dopar% {
    single_comparison(sample_size = sample_size, i = i)
  }
  
  stopCluster(cluster)
  print("Simulation complete.")
  return(as.data.frame(simulation_data))
}

# Number of iterations
nsim <- 200

# Sample sizes - list
sample_sizes <- c(sample_size/2, sample_size, 3*sample_size/2)

res <- NULL
for (n in sample_sizes) {
  res <- rbind(res, cbind(simulation(n, nsim), n))
}

metrics <- c("cs", "auc", "brier", "mape")
metric_labels <- c(cs = "Calibration Slope", auc = "AUC", brier = "Brier Score", mape = "MAPE")

for (m in metrics) {
  # Scatter Plot
  res_m_long <- res %>%
    select(n, starts_with(m)) %>%
    rename_with(~str_remove(., paste0(m, "_")), starts_with(m)) %>%
    pivot_longer(cols = c(boot, split, cv), names_to = "method", values_to = "internal")
    
  cors_m <- res_m_long %>%
    group_by(n, method) %>%
    summarise(cor_val = cor(internal, true, method = "spearman"), .groups = "drop")
    
  scatter <- ggplot(res_m_long, aes(x = internal, y = true, color = method)) +
    geom_point(alpha = 0.6) +
    facet_wrap(~ n + method, scales = "free") +
    theme_minimal() +
    geom_text(
      data = cors_m,
      aes(x = -Inf, y = Inf, label = paste0("r = ", round(cor_val, 3))),
      hjust = -0.1, vjust = 1.2, inherit.aes = FALSE, size = 4
    ) +
    labs(
      title = paste0(metric_labels[m], " correlation (internal vs external 'true')"),
      subtitle = paste0("Sims = ", nsim),
      x = paste0("Internal ", metric_labels[m], " estimate"),
      y = paste0("External 'true' ", metric_labels[m])
    )
    
  ggsave(paste0("scatter_", m, ".png"), scatter, width = 10, height = 8)
  
  # Box Plot
  res_m_box <- res %>%
    select(n, starts_with(m)) %>%
    rename_with(~str_remove(., paste0(m, "_")), starts_with(m)) %>%
    pivot_longer(cols = -n, names_to = "method", values_to = "val")
    
  box <- ggplot(res_m_box, aes(x = method, y = val, fill = method)) +
    geom_boxplot(alpha = 0.7) +
    theme_minimal() + 
    facet_wrap(~n) +
    labs(
      title = paste0("Comparison of ", metric_labels[m]),
      x = "",
      y = metric_labels[m]
    ) +
    theme(legend.position = "none")
    
  # Add horizontal line for reference if applicable
  if (m == "cs") box <- box + geom_hline(yintercept = 0.9, color = "blue", linetype = "dotted")
  if (m == "auc") box <- box + geom_hline(yintercept = auc_val, color = "blue", linetype = "dotted")
  
  ggsave(paste0("box_", m, ".png"), box, width = 10, height = 8)
}

print("All plots saved.")