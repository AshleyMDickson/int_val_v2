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

# True model performance - Cal Slope only for now
# test_data must be supplied as a dataframe; note that dgp() returns a list with df and auc_true, so we need to use test_data$df when calling this function
model_truth <- function(df, test_data) {
  model <- glm(outcome ~ ., family = "binomial", data = df)
  ext_lp <- predict(model, newdata = test_data, type = "link") # test_data drawn once outside this function
  cs_true <- unname(coef(glm(test_data$outcome ~ ext_lp, family = "binomial"))[2])
  #auc_true <- suppressMessages(auc(response = test_data$df$outcome, predictor = ext_lp))
  #brier_true <- mean((ext_lp - test_data$df$outcome)^2)
  #mape_true <- mean(abs(ext_lp - test_data$df$outcome))
  return(cs_true)
}
## Trial run to make sure it works
# dev_data <- dgp(sample_size, k, return_auc = FALSE)
# apparent_slope <- model_truth(df=dev_data, test_data = dev_data)
# print(paste0("Apparent calibration slope = ", apparent_slope))
# print(paste0("True calibration slope = ", model_truth(df=dev_data, test_data = test_data$df)))

# Get Bootstrap corrected slope = apparent[= 1.0] - optimism[= 1.0 - mean_test_slope]
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
#print(boot_corr_slope(dev_data, B = 200))

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
#print(sample_split_slope(dev_data, prop_train = 0.5))
#x <- numeric(1000)
#for (i in 1:1000) {x[i] <- sample_split_slope(dev_data, prop_train = 0.8)}
#hist(x, breaks = 20)

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
#y <- numeric(1000)
#for (i in 1:1000) {y[i] <- cv_slope(dev_data, K = 10)}
#hist(y, breaks = 20)

# See MP's code for "cv_slope_try" for an oversampled CV slope calculation. 
# This enables model to be trained on a dataset of full size.

single_comparison <- function(sample_size, i) {
  set.seed(i)
  df <- dgp(sample_size, k)
  test_df <- dgp(100000, k)
  
  cs_true  <- model_truth(df, test_df)
  cs_boot  <- boot_corr_slope(df, B = 200)
  cs_split <- sample_split_slope(df, prop_train = 0.5)
  cs_cv    <- cv_slope(df, K = 10)

  c(cs_true = cs_true, cs_boot = cs_boot, cs_split = cs_split, cs_cv = cs_cv)
}
#print(single_comparison(sample_size, 1))

simulation <- function(sample_size, nsim) { # n = number of simualtion repetitions
  cores <- detectCores() - 1
  cluster <- makeCluster(cores)
  registerDoParallel(cluster)
  print(paste0("Running simulation on ", cores, " cores..."))
  
  simulation_data <- foreach(
    i = 1:nsim,
    .combine = rbind,
    .packages = c("rms", "pROC"),
    .export = c("single_comparison", "dgp", "model_truth", "cv_slope", 
    "sample_split_slope",
                "boot_corr_slope", "alpha", "beta", "k")
  ) %dopar% {
    single_comparison(sample_size = sample_size, i = i)
  }
  
  stopCluster(cluster)
  #simulation_data <- t(replicate(n, single_comparison()))
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
#res


cors <- res %>%
  group_by(n) %>%
  summarise(
    cor_val1 = cor(cs_boot, cs_true, method = "spearman"),
    cor_val2 = cor(cs_cv, cs_true, method = "spearman"),
    cor_val3 = cor(cs_split, cs_true, method = "spearman")
      )

cors_long <- cors %>%
  pivot_longer(
    cols = starts_with("cor_val"),
    names_to = "method",
    values_to = "cor_val"
  ) %>%
  mutate(
    # Match the column names in res_long
    method = case_when(
      method == "cor_val1" ~ "cs_boot",
      method == "cor_val2" ~ "cs_cv",
      method == "cor_val3" ~ "cs_split"
    )
  )

res_long <- res %>%
  pivot_longer(cols = c(cs_boot, cs_split, cs_cv),
               names_to = "method", values_to = "cs_internal")


scatter <- ggplot(res_long, aes(x = cs_internal, y = cs_true, color = method)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ n + method, scales = "free_x") +
  theme_minimal() +
  coord_cartesian(xlim = c(0.4, 1.6),
                  ylim = c(0.4, 1.6)) +
  geom_text(
    data = cors_long,
    aes(
      x = -Inf,
      y = Inf,
      label = paste0("r = ", round(cor_val, 3))
    ),
    hjust = -0.1,
    vjust = 1.2,
    inherit.aes = FALSE,
    size = 4
  )+
  labs(
    title = "Calibration slope correlation (internal vs external 'true')",
    subtitle = paste0(
      "Sims = ", nsim
    ),
    x = "Internal calibration slope estimate",
    y = "External 'true' calibration slope"
  )

scatter

ggsave("scatter.png")



# Boxplot 
# Reshape to long format
res_box <- res %>%
  pivot_longer(
    cols = c(cs_true, cs_boot, cs_cv, cs_split),
    names_to = "method",
    values_to = "cal_slope"
  )

box <- ggplot(res_box, aes(x = method, y = cal_slope, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() + facet_wrap(~n)+
  labs(
    title = "Comparison of Calibration Slopes",
    x = "",
    y = "Calibration Slope"
  ) +
  theme(legend.position = "none") +  # Horizontal mean line
  geom_hline(
    aes(yintercept =0.9),
    color = "blue",
    linetype = "dotted",
    linewidth = 1
  ) 

box

ggsave("box.png")

# cors <- c(
#   boot  = cor(res$cs_boot,  res$cs_true),
#   split = cor(res$cs_split, res$cs_true),
#   cv    = cor(res$cs_cv,    res$cs_true)
# )
# print(cors)

# res_long <- res %>%
#   pivot_longer(cols = c(cs_boot, cs_split, cs_cv),
#                names_to = "method", values_to = "cs_internal")

# plt <- ggplot(res_long, aes(x = cs_internal, y = cs_true)) +
#   geom_jitter(width = 0.5, height = 0.5, alpha = 0.5) +
#   facet_wrap(~ method, scales = "free_x") +
#   theme_minimal() +
#   labs(
#     title = "Calibration slope correlation (internal vs external 'true')",
#     subtitle = paste0(
#       "Sample size = ", sample_size,
#       ". Sims = ", nsim,
#       ". cor(boot)=", round(cors["boot"], 3),
#       ", cor(split)=", round(cors["split"], 3),
#       ", cor(cv)=", round(cors["cv"], 3)
#     ),
#     x = "Internal calibration slope estimate",
#     y = "External 'true' calibration slope"
#   )

# print(plt)