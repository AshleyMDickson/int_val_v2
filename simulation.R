# Import libraries
library(pROC)
library(rms)
library(samplesizedev)
library(foreach)
library(doParallel)

# # Set-up parameters and such
#sampsizedev <- suppressWarnings(samplesizedev(outcome = "Binary", S = 0.9, phi = 0.15, c = 0.85, p = 10))
#sample_size <- sampsizedev$sim # tried up to 10,240 sample size
sample_size <- 631
beta <- c(0.45, 0.40, -0.35, 0.30, -0.25, 1.20, 0.15, 0.10, -1.8, 0.05)
#set.seed(1729)

get_alpha <- function(target_prev, beta, k) {
    N <- 100000
    predictors <- matrix(rnorm(N*k), nrow = N, ncol = k)
    lp_fixed <- predictors %*% beta
    root <- function(a) mean(plogis(a + lp_fixed)) - target_prev
    alpha <- uniroot(root, c(-10, 10))$root
    return(alpha)
}

alpha <- get_alpha(0.15, beta, length(beta))

dgp <- function(n, k) {
    predictors <- matrix(rnorm(n*k), nrow = n, ncol = k)
    linear_predictor <- alpha + (predictors %*% beta)
    probability <- 1/(1+exp(-1*linear_predictor))
    outcome <- rbinom(n, size = 1, prob = probability)
    auc(outcome, probability)
    return(data.frame(outcome, predictors))
}

# true model performance - Cal Slope only for now
test_data <- dgp(200000, length(beta))

model_truth <- function(df, test_data) {
    #df <- dgp(sample_size, length(beta))
    model <- glm(outcome ~ ., family = "binomial", data = df)
    ## check apparent slope = 1 to ensure correct calculation
    #pred_lp <- predict(model, newdata = df, type = "link")
    #cs_naive <- unname(coef(glm(df$outcome ~ pred_lp, family = "binomial"))[2])

    #new_data <- dgp(200000, length(beta)) # created externally now for speed
    ext_lp <- predict(model, newdata = test_data, "link")
    cs_true <- unname(coef(glm(test_data$outcome ~ ext_lp, family = "binomial"))[2])
    
    #auc_true <- suppressMessages(auc(response = new_data$outcome, predictor = ext_preds))
    #brier_true <- mean((ext_preds - new_data$outcome)^2)
    #mape_true <- mean(abs(ext_preds - new_data$outcome))

    return(c(cs_true = cs_true))
}
#res <- model_truth()
#print(res)

# get Bootstrap corrected slope = apparent[= 1.0] - optimism[= 1.0 - mean_test_slope]
# ? wrong ?
boot_corr_slope <- function(df, B = 200) {
    bootstrap_slopes <- numeric(B)
    n <- nrow(df)

    for(i in 1:B) {
        index_resampled <- sample(n, replace = TRUE)
        bootstrap_df <- df[index_resampled,] # resampled df

        bootstrap_model <- glm(outcome ~ ., family = "binomial", data = bootstrap_df)
        predictions_orig <- predict(bootstrap_model, newdata = df, type = "link") # get Bootstrap model preds on original data
        bootstrap_slopes[i] <- coef(glm(df$outcome ~ predictions_orig, family = "binomial"))[2] # individual Bootstrap cal. slope
    }
    return(mean(bootstrap_slopes))
}

single_comparison <- function() {
    df <- dgp(sample_size, length(beta)) # training data
    # model <- glm(outcome ~ ., family = "binomial", data = df) # main model
    
    # ext_data <- dgp(100000, length(beta)) # large test dataset
    # ext_lp <- predict(model, newdata = ext_data, type = "link")
    # cs_true <- unname(coef(glm(ext_data$outcome ~ ext_lp, family = "binomial"))[2]) 
    cs_true <- model_truth(df, test_data) # true cal slope
    boot <- boot_corr_slope_v2(df = df, B = 200) # optimism-corrected slope

    cs_app <- boot$apparent_slope_original
    opt_boot <- boot$mean_optimism
    cs_bootcorr <- boot$corrected_slope

    opt_true <- cs_app - cs_true 

    return(
        c(
            True = cs_true, 
            App = cs_app,
            Boot = cs_bootcorr,
            OptTrue = opt_true,
            OptBoot = opt_boot
        )
    )
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
        "boot_corr_slope_v2", "alpha", "beta", "sample_size", "test_data")) %dopar% {
            single_comparison()
        }

    stopCluster(cluster)
    print("Simulation complete.")
    return(as.data.frame(simulation_data))
}

# res <- simulation(500)
# str(res)
# summary(res)

summary_stats <- c(
    n = nrow(res),

    mean_true = mean(res$True),
    sd_true   = sd(res$True),

    mean_app  = mean(res$App),
    sd_app    = sd(res$App),

    mean_boot = mean(res$Boot),
    sd_boot   = sd(res$Boot),

    cor_true_boot = cor(res$True, res$Boot),

    mean_opt_true = mean(res$OptTrue),
    sd_opt_true   = sd(res$OptTrue),

    mean_opt_boot = mean(res$OptBoot),
    sd_opt_boot   = sd(res$OptBoot),

    cor_opt_true_opt_boot = cor(res$OptTrue, res$OptBoot),

    range = range(res$Boot),
    unique_vals = length(unique(round(res$Boot, 3)))
)

#print(summary_stats)

# png("diag1.png")
# plot(res$True, res$Boot, pch=16, col=rgb(0,0,0,0.3))
# abline(0,1,lwd=2)
# dev.off()

# png("diag2.png")
# plot(res$OptTrue, res$OptBoot, pch=16, col=rgb(0,0,0,0.3))
# abline(0,1,lwd=2)
# dev.off()

