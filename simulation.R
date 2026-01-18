# Import libraries
library(pROC)
library(rms)
library(samplesizedev)
library(foreach)
library(doParallel)

# # Set-up parameters and such
# sampsizedev <- suppressWarnings(samplesizedev(outcome = "Binary", S = 0.9, phi = 0.15, c = 0.85, p = 10))
# sample_size <- sampsizedev$sim # tried up to 10,240 sample size
sample_size <- 631
beta <- c(0.45, 0.40, -0.35, 0.30, -0.25, 1.20, 0.15, 0.10, -1.8, 0.05)
set.seed(1729)

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

truth_convergence_plot <- function(n) {

    set.seed(1729)
    n_sims <- n
    pb <- txtProgressBar(min = 0, max = n_sims, style = 3) # progress bar
    raw_results <- numeric(n_sims) # loop it
    for(i in 1:n_sims) {
        df <- dgp(sample_size, length(beta))
        raw_results[i] <- model_truth(df)
        setTxtProgressBar(pb, i)
    }
    close(pb)

    sim_results <- data.frame(cs_true = raw_results)
    sequence_n <- seq_along(sim_results$cs_true)

    sim_results$cum_mean <- cumsum(sim_results$cs_true) / sequence_n
    cum_sq <- cumsum(sim_results$cs_true^2) # for cum_var
    cum_var <- (cum_sq - sequence_n*sim_results$cum_mean^2)/(sequence_n-1)
    cum_var[1] <- 0 # 1st element set to 0 as no var

    #sim_results$cum_sd <- sapply(sequence_n, function(x) {
    #    if (x == 1) return(0)
    #    sd(sim_results$cs_true[1:x])
    #}) #removed as inefficient

    sim_results$cum_sd <- sqrt(cum_var)
    sim_results$cum_se <- sim_results$cum_sd/sqrt(sequence_n)
    sim_results$lower_ci <- sim_results$cum_mean - 1.96*sim_results$cum_se
    sim_results$upper_ci <- sim_results$cum_mean + 1.96*sim_results$cum_se

    png("calibration_convergence.png", width = 2000, height = 1500, res = 300)

    y_lims <- c(0.8, 1.05)

    plot(x = 1:n_sims, y = sim_results$cum_mean, type = "n", ylim = y_lims, # type n = draw nothing yet
        xlab = "No. of simulations", ylab = "Cumulative mean of True Slope",
        main = "Convergence to Expected Calibration Slope",
        las = 1)

    polygon(
        x = c(sequence_n, rev(sequence_n)),
        y = c(sim_results$lower_ci, rev(sim_results$upper_ci)),
        col = rgb(0.8, 0.8, 0.8, 0.5),
        border = NA
    )
    lines(x = sequence_n, y = sim_results$cum_mean, lwd = 2)
    abline(h = tail(sim_results$cum_mean, 1), col = "red", lty = 2)
    dev.off()
    print("Convergence plot saved.")
}
#truth_convergence_plot(2000)

# get Bootstrap corrected slope = apparent[= 1.0] - optimism[= 1.0 - mean_test_slope]
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
    cs_boot <- boot_corr_slope(df = df, B = 200) # optimism-corrected slope
    return(c(True = cs_true, Boot = cs_boot))
}

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
        "boot_corr_slope", "alpha", "beta", "sample_size", "test_data")) %dopar% {
            single_comparison()
        }

    stopCluster(cluster)
    print("Simulation complete.")
    return(as.data.frame(simulation_data))
}

boxplot_correlation <- function() {
    sim_results <- simulation(200)
    png("calibration_comparison.png", width = 2400, height = 1200, res = 300)
    par(mfrow = c(1, 2), mar = c(5, 5, 4, 2)) # 1 row, 2 cols

    boxplot( #plot 1 - LHS
        sim_results,
        main = "Comparison of Model Calibration Slopes \n(True Slope vs Bootstrap Estimate)",
        ylab = "Calibration Slope",
        col = c("#69b3a2", "#404080"),
        ylim = c(0.6, 1.4),
        las = 1
        )

    abline(h = 1.0, col = "grey", lty = 3) # apparent slope reference line
    points(1:2, colMeans(sim_results), pch = 18, col = "#f3a1a1", cex = 2) # mean markers

    plot( # plot 2 - RHS
        x = sim_results$True, y = sim_results$Boot, # check: is this the right way round?
        pch = 16, col = rgb(0, 0, 0, 0.4),
        xlim = c(0.6, 1.4), ylim = c(0.6, 1.4),
        xlab = "True Slope (external)", ylab = "Bootstrap Est. Slope",
        main = "Individual Agreement\n(1 point = 1 simulation run)",
        las = 1)

    abline(0, 1, col = "#a8f0e3", lwd = 2) # perfect agreement
    abline(h = 1, v = 1, col = "grey", lty = 2)

    r_val <- cor(sim_results$True, sim_results$Boot)
    text(0.7, 1.3, paste0("Correlation = ", round(r_val, 2)), pos = 4, col = "grey")

    dev.off()
    print("Boxplot and correlation saved.")
    }
#boxplot_correlation()

correlation_dist <- function(n) {
    r_vals <- numeric(n)
    for (i in 1:n) {
        res <- simulation(n)
        r_vals[i] <- cor(res$True, res$Boot)
        print(paste0("Finished batch ", i))
    }

    png("correlation_distribution.png")
    hist(
        r_vals,
        main = "Distribution of True-Boot Correlations", xlab = "Correlation Coefficient"
        width = 2400, height = 1800, res = 300
    )
    dev.off()
}
correlation_dist(100)