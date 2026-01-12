# Import libraries
library(pROC)
library(rms)
library(samplesizedev)

# # Set-up parameters and such
# sampsizedev <- suppressWarnings(samplesizedev(outcome = "Binary", S = 0.9, phi = 0.15, c = 0.85, p = 10))
# sample_size <- sampsizedev$sim
sample_szie <- 631

beta <- c(0.45, 0.40, -0.35, 0.30, -0.25, 1.20, 0.15, 0.10, -1.8, 0.05)
set.seed(1729)

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

# true model performance
model_truth <- function() {
    df <- dgp(sample_size, length(beta))
    model <- glm(outcome ~ ., family = "binomial", data = df)
    new_data <- dgp(1000000, length(beta))
    ext_preds <- predict(model, new_data, "response")
    ext_lp <- qlogis(ext_preds)

    auc_true <- suppressMessages(auc(response = new_data$outcome, predictor = ext_preds))
    cs_true <- coef(glm(new_data$outcome ~ ext_lp, family = "binomial"))[2]
    brier_true <- mean((ext_preds - new_data$outcome)^2)
    mape_true <- mean(abs(ext_preds - new_data$outcome))
    return(list(auc_true, cs_true, brier_true, mape_true))
}

res <- model_truth()
print(res)