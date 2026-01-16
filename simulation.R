# Import libraries
library(pROC)
library(rms)
library(samplesizedev)

# # Set-up parameters and such
# sampsizedev <- suppressWarnings(samplesizedev(outcome = "Binary", S = 0.9, phi = 0.15, c = 0.85, p = 10))
# sample_size <- sampsizedev$sim
sample_size <- 631

beta <- c(0.45, 0.40, -0.35, 0.30, -0.25, 1.20, 0.15, 0.10, -1.8, 0.05)
#set.seed(1729)

get_alpha <- function(target_prev, beta, k) {
    N <- 1000000
    predictors <- matrix(rnorm(N*k), nrow = N, ncol = k)
    lp_fixed <- predictors %*% beta
    root <- function(a) mean(plogis(a + lp_fixed)) - target_prev
    alpha <- uniroot(root, c(-10, 10))$root
    return(alpha)
}

alpha <- get_alpha(0.15, beta, 10)

dgp <- function(n, k) {
    predictors <- replicate(k, rnorm(n))
    linear_predictor <- alpha + (predictors %*% beta)
    probability <- 1/(1+exp(-1*linear_predictor))
    outcome <- rbinom(n, size = 1, prob = probability)
    return(data.frame(outcome, predictors))
}

# true model performance - Cal Slope only for now
model_truth <- function() {
    df <- dgp(sample_size, length(beta))
    model <- glm(outcome ~ ., family = "binomial", data = df)
    ## check naive slope = 1 to ensure correct calculation
    #pred_lp <- predict(model, newdata = df, type = "link")
    #cs_naive <- unname(coef(glm(df$outcome ~ pred_lp, family = "binomial"))[2])

    new_data <- dgp(200000, length(beta))
    ext_lp <- predict(model, newdata = new_data, "link")
    cs_true <- unname(coef(glm(new_data$outcome ~ ext_lp, family = "binomial"))[2])
    
    #auc_true <- suppressMessages(auc(response = new_data$outcome, predictor = ext_preds))
    #brier_true <- mean((ext_preds - new_data$outcome)^2)
    #mape_true <- mean(abs(ext_preds - new_data$outcome))

    return(c(cs_true = cs_true))
}

#res <- model_truth()
#print(res)

set.seed(1729)
n_sims <- 100
raw_results <- replicate(n_sims, model_truth())
sim_results <- data.frame(cs_true = raw_results)

sim_results$cum_mean <- cumsum(sim_results$cs_true) / seq_along(sim_results$cs_true)

png("calibration_convergence.png", width = 2000, height = 1500, res = 300)

plot(x = 1:n_sims, y = sim_results$cum_mean, type = "l", ylim = c(0.8, 1.0),
    xlab = "No. of simulations", ylab = "Cumulative mean of True Slope",
    main = "Convergence to Expected Calibration Slope",
    las = 1)
abline(h = tail(sim_results$cum_mean, 1), col = "red", lty = 2)

dev.off()