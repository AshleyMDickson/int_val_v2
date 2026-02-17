# Import libraries
library(pROC)
#library(rms)
library(samplesizedev)
library(foreach)
library(doParallel)
library(ggplot2)
library(tidyverse)

# # Set-up parameters and such
beta <- c(0.45, 0.40, -0.35, 0.30, -0.25, 1.20, 0.15, 0.10, -1.8, 0.05)
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

#alpha <- get_alpha(phi, beta, k)
alpha <- -3

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
# test_data <- dgp(100000, k, return_auc = TRUE)
# auc_val <- test_data$auc_true
#samp_size_calc <- suppressWarnings(samplesizedev(outcome = "Binary", S = 0.9, phi = phi, c = auc_val, p = k))
#sample_size <- samp_size_calc$sim
sample_size <- 592 #592
sims <- 200

# true model performance - Cal Slope only for now
model_truth <- function(df, test_data) {
  model <- glm(outcome ~ ., family = "binomial", data = df)
  ext_lp <- predict(model, newdata = test_data, type = "link") # test_data drawn once outside this function
  cs_true <- unname(coef(glm(test_data$outcome ~ ext_lp, family = "binomial"))[2])
  #auc_true <- suppressMessages(auc(response = new_data$outcome, predictor = ext_preds))
  #brier_true <- mean((ext_preds - new_data$outcome)^2)
  #mape_true <- mean(abs(ext_preds - new_data$outcome))
  return(cs_true)
}

# get Bootstrap corrected slope = apparent[= 1.0] - optimism[= 1.0 - mean_test_slope]
# Here, the goal is to try a couple of things: (1) try getting predictions not in original data but in another bootstrap resample, and (2) double bootstrap.
# What does double bootstrap mean?
# Oh and also, Ambler's empirical-Bayes thing. What does it mean?!?!?

boot_corr_slope <- function(df, B = 200) {
  bootstrap_slopes <- numeric(B)
  n <- nrow(df)
  
  for(i in 1:B) {
    index_resampled <- sample(n, replace = TRUE)
    bootstrap_df <- df[index_resampled,] # resampled df
    bootstrap_model <- glm(outcome ~ ., family = "binomial", data = bootstrap_df)
    # what about test in orig data less bootstrap sampled records? is this 0.632????
    #oob_indices <- which(!(1:n %in% index_resampled))
    #df_not_bootstrap <- df[sample(oob_indices, size = n, replace = TRUE), ] # original df, less resampled records
    #predictions_not_boot <- predict(bootstrap_model, newdata = df_not_bootstrap, type = "link")
    predictions_orig <- predict(bootstrap_model, newdata = df, type = "link") # get Bootstrap model preds on original data
    bootstrap_slopes[i] <- coef(glm(df$outcome ~ predictions_orig, family = "binomial"))[2] # individual Bootstrap cal. slope
  }
  return(mean(bootstrap_slopes)) 
  # what about trying a random sample of the 200 bootstrap slopes? or max? or min? median? something NOT mean!
  #return(sample(bootstrap_slopes, 1)) # random sample of bootstrap slopes
}
#boot_cs <- boot_corr_slope(B = 200)
#print(paste0("Bootstrap-corrected Calibration Slope = ", boot_cs))

sample_split_slope <- function(df, prop_train = 0.5) {
  n <- nrow(df)
  idx <- sample.int(n)
  n_train <- floor(prop_train * n)
  train_idx <- idx[1:n_train]
  test_idx  <- idx[(n_train + 1):n]
  
  train_df <- df[train_idx, , drop = FALSE]
  test_df  <- df[test_idx,  , drop = FALSE]
  
  model <- glm(outcome ~ ., family = "binomial", data = train_df)
  lp_test <- predict(model, newdata = test_df, type = "link")
  unname(coef(glm(test_df$outcome ~ lp_test, family = "binomial"))[2])
}

cv_slope <- function(df, K = 10) {
  n <- nrow(df)
  folds <- sample(rep(1:K, length.out = n))
  lp_oof <- rep(NA_real_, n) # pooled out of fold linear preds
  
  for (kfold in 1:K) {
    train_df <- df[folds != kfold, , drop = FALSE]
    test_df  <- df[folds == kfold, , drop = FALSE]
    model <- glm(outcome ~ ., family = "binomial", data = train_df)
    lp_oof[folds == kfold] <- predict(model, newdata = test_df, type = "link")
  }
  unname(coef(glm(df$outcome ~ lp_oof, family = "binomial"))[2])
}

single_comparison <- function() {
  df <- dgp(sample_size, k)
  test_df <- dgp(100000, k)
  
  cs_true  <- model_truth(df, test_df)
  cs_boot  <- boot_corr_slope(df, B = 200)
  cs_split <- sample_split_slope(df, prop_train = 0.5)
  cs_cv    <- cv_slope(df, K = 10)
  
  c(cs_true = cs_true, cs_boot = cs_boot, cs_split = cs_split, cs_cv = cs_cv)
}
#print(single_comparison())

simulation <- function(n) { # n = number of simualtion repetitions
  cores <- detectCores() - 1
  cluster <- makeCluster(cores)
  registerDoParallel(cluster)
  print(paste0("Running simulation on ", cores, " cores..."))
  
  simulation_data <- foreach(
    i = 1:n,
    .combine = rbind,
    .packages = c("rms", "pROC"),
    .export = c("single_comparison", "dgp", "model_truth", "cv_slope", "sample_split_slope",
                "boot_corr_slope", "alpha", "beta", "k","sample_size")) %dopar% {
                  single_comparison()
                }
  
  stopCluster(cluster)
  #simulation_data <- t(replicate(n, single_comparison()))
  print("Simulation complete.")
  return(as.data.frame(simulation_data))
}

res <- simulation(sims)

cors <- c(
  boot  = cor(res$cs_boot,  res$cs_true),
  split = cor(res$cs_split, res$cs_true),
  cv    = cor(res$cs_cv,    res$cs_true)
)
print(cors)

res_long <- res %>%
  pivot_longer(cols = c(cs_boot, cs_split, cs_cv),
               names_to = "method", values_to = "cs_internal")

plt <- ggplot(res_long, aes(x = cs_internal, y = cs_true)) +
  geom_jitter(width = 0.5, height = 0.5, alpha = 0.5) +
  facet_wrap(~ method, scales = "free_x") +
  theme_minimal() +
  labs(
    title = "Calibration slope correlation (internal vs external 'true')",
    subtitle = paste0(
      "Sample size = ", sample_size,
      ". Sims = ", sims,
      ". cor(boot)=", round(cors["boot"], 3),
      ", cor(split)=", round(cors["split"], 3),
      ", cor(cv)=", round(cors["cv"], 3)
    ),
    x = "Internal calibration slope estimate",
    y = "External 'true' calibration slope"
  )

print(plt)