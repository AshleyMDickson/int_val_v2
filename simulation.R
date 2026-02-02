# Import libraries
library(pROC)
#library(rms)
#library(samplesizedev)
library(foreach)
library(doParallel)

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

alpha <- get_alpha(phi, beta, k)

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

test_data <- dgp(100000, k, return_auc = TRUE)
test_df <- test_data$df
auc_val <- test_data$auc_true
#print(auc_val)

#sampsizedev <- suppressWarnings(samplesizedev(outcome = "Binary", S = 0.9, phi = phi, c = auc_val, p = k))
#sample_size <- sampsizedev$sim
sample_size <- 592

# true model performance - Cal Slope only for now
model_truth <- function(test_data) {
    df <- dgp(sample_size, k) # Model data with sample size, re-drawn inside function
    model <- glm(outcome ~ ., family = "binomial", data = df)
    ext_lp <- predict(model, newdata = test_data, type = "link") # test_data drawn once outside this function
    cs_true <- unname(coef(glm(test_data$outcome ~ ext_lp, family = "binomial"))[2])
    
    #auc_true <- suppressMessages(auc(response = new_data$outcome, predictor = ext_preds))
    #brier_true <- mean((ext_preds - new_data$outcome)^2)
    #mape_true <- mean(abs(ext_preds - new_data$outcome))
    return(cs_true)
}

#true_cs <- model_truth(test_df)
#print(paste0("True Calibration Slope = ", true_cs))
#hist(true_cs)
#print(paste0("Mean=",mean(true_cs), "; SD=", sd(true_cs)))

# get Bootstrap corrected slope = apparent[= 1.0] - optimism[= 1.0 - mean_test_slope]
# Here, the goal is to try a couple of things: (1) try getting predictions not in original data but in another bootstrap resample, and (2) double bootstrap.
# What does double bootstrap mean?

# Oh and also, Ambler's empirical-Bayes thing. What does it mean?!?!?

boot_corr_slope <- function(B = 200) {
    bootstrap_slopes <- numeric(B)
    df <- dgp(sample_size, k) # Model data with sample size, re-drawn inside function
    n <- nrow(df)

    for(i in 1:B) {
        index_resampled <- sample(n, replace = TRUE)
        bootstrap_df <- df[index_resampled,] # resampled df
        bootstrap_model <- glm(outcome ~ ., family = "binomial", data = bootstrap_df)

        #index_resample_2
        predictions_orig <- predict(bootstrap_model, newdata = df, type = "link") # get Bootstrap model preds on original data
        bootstrap_slopes[i] <- coef(glm(df$outcome ~ predictions_orig, family = "binomial"))[2] # individual Bootstrap cal. slope
    }
    # return(mean(bootstrap_slopes)) # what about trying a random sample of the 200 bootstrap slopes? or max? or min? median? something NOT mean!
    return(max(bootstrap_slopes))
}
#boot_cs <- boot_corr_slope(B = 200)
#print(paste0("Bootstrap-corrected Calibration Slope = ", boot_cs))

single_comparison <- function() {
    cs_true <- model_truth(test_df) # true cal slope
    cs_boot <- boot_corr_slope(B = 200) # optimism-corrected slope
    return(c(cs_true = cs_true, cs_boot = cs_boot))
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
        .export = c("single_comparison", "dgp", "model_truth",
        "boot_corr_slope", "alpha", "beta", "k","sample_size", "test_df")) %dopar% {
            single_comparison()
        }

    stopCluster(cluster)
    #simulation_data <- t(replicate(n, single_comparison()))
    print("Simulation complete.")
    return(as.data.frame(simulation_data))
}

res <- simulation(100)
plot(res)
