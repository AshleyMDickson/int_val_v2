# not really needed now
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

# not needed in the end
boot_corr_slope_v2 <- function(df, B = 200) {
    model_orig <- glm(outcome ~ ., family = "binomial", data = df) # original model
    lp_orig <- predict(model_orig, newdata = df, type = "link")
    app_cs_orig <- unname(coef(glm(df$outcome ~ lp_orig, family = "binomial"))[2]) # apparent slope

    n <- nrow(df)

    app_orig_vec <- rep(app_cs_orig, B) # same value, B repeats
    app_boot_vec <- numeric(B)
    test_boot_vec <- numeric(B)
    optimism_vec <- numeric(B)

    for (i in 1:B) {
        index_resampled <- sample(n, replace = TRUE)
        bootstrap_df <- df[index_resampled,] # resamnpled df
        model_bootstrap <- glm(outcome ~ ., family = "binomial", data = bootstrap_df)
        lp_boot_in_boot <- predict(model_bootstrap, newdata = bootstrap_df, type = "link")
        app_cs_bootstrap <- unname(coef(glm(outcome ~ lp_boot_in_boot, family = "binomial", data = bootstrap_df))[2]) # apparenbt bootstrap cs

        lp_boot_in_orig <- predict(model_bootstrap, newdata = df, type = "link")
        test_cs_boot <- unname(coef(glm(outcome ~ lp_boot_in_orig, family = "binomial", data = df))[2]) # model cal slope on bootstrap data

        app_boot_vec[i] <- app_cs_bootstrap
        test_boot_vec[i] <- test_cs_boot
        optimism_vec[i] <- app_cs_bootstrap - test_cs_boot
    }

    steps <- data.frame(
        app_slope_orig = app_orig_vec,
        app_slope_boot = app_boot_vec,
        test_slope_boot = test_boot_vec
    )

    cs_corrected <- app_cs_orig - mean(optimism_vec)
    
    list(
        corrected_slope = cs_corrected, # single value
        apparent_slope_original = app_cs_orig,
        mean_optimism = mean(optimism_vec),
        steps = steps,
        optimism = optimism_vec
    )
}

# boxplot and correlation figure
boxplot_correlation <- function() {
    sim_results <- simulation(200)
    png("calibration_comparison_v2.png", width = 2400, height = 1200, res = 300)
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

# histogram of correlations
correlation_dist <- function(n) {
    r_vals <- numeric(n)
    for (i in 1:n) {
        res <- simulation(n)
        r_vals[i] <- cor(res$True, res$Boot)
        print(paste0("Finished batch ", i))
    }

    png("correlation_distribution_v2.png", width = 2400, height = 1800, res = 300)
    hist(
        r_vals,
        main = "Distribution of True-Boot Correlations", xlab = "Correlation Coefficient"
    )
    dev.off()
    
}
#correlation_dist(10)