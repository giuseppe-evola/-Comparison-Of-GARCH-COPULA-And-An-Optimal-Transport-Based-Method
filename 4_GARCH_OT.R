# ==============================================================================
#                               GARCH-OT Model
# ==============================================================================
# This script implements the algorithm used to compute the one-step-ahead 
# Value at Risk (VaR at T+1) based on a GARCH-OT model with static dependence 
# structure (i.e., the OT model is fitted only once on the training set).
#
# Before that, we provide the charts to further validate the adherence of the simulated data to the real ones.
#
# The procedure operates as follows:
#
# 1. Define rolling windows of 1500 observations for each time series using 
#    the `rollapply` function.
#
# 2. Fit the OT model to the training data. 
#    → This step is performed externally in MATLAB (see folder: 'OT ALGORITHM definitive').
#    → Here, we simply import precomputed datasets containing 1,000,000×4 standardized 
#      residual simulations for each specification. 
#    → Three datasets have been imported, corresponding to the top three OT specifications.
#
# 3. For each window, apply the best marginal model for each stock index, 
#    specifically ARMA(0,0) - TGARCH(1,1) with skewed t-distributed innovations.
#
# 4. Forecast the one-step-ahead conditional volatility for each index.
#
# 5. Transform the standardized residuals into standard uniform variables 
#    (pseudo-observations) using the Probability Integral Transform (PIT).
#
# 6. For each OT specification, sample 10,000x4 pseudo-observations.
#
# 7. Transform the simulated pseudo-observations into simulated log-returns by 
#    multiplying them by the forecasted volatilities.
#
# 8. Aggregate the simulated log-returns to obtain portfolio log-returns and 
#    convert them into arithmetic returns.
#
# 9. Estimate the Value at Risk by computing the 5th and 1st percentiles of the 
#    simulated distribution of portfolio returns.
#
# Steps 3–9 are repeated for each rolling window over the test set, which consists 
# of 1,110 observations.
#
# The entire algorithm is parallelized to improve computational efficiency.
# ==============================================================================

setwd("C:/Users/giuse/OneDrive/Desktop/Preparazione tesi/Thesis Code")


# Import the dataframes of the three best specifications (epsilon-bins)
sinkhorn_simulations_TOP_1 <- read.csv("OT Algorithm definitive/sinkhorn_simulations_TOP1.csv")
sinkhorn_simulations_TOP_2 <- read.csv("OT Algorithm definitive/sinkhorn_simulations_TOP2.csv")
sinkhorn_simulations_TOP_3 <- read.csv("OT Algorithm definitive/sinkhorn_simulations_TOP3.csv")

sim_list <- list(
  OT_model_1 = sinkhorn_simulations_TOP_1,  
  OT_model_2 = sinkhorn_simulations_TOP_2,
  OT_model_3 = sinkhorn_simulations_TOP_3)

param_tuning_summary <- read.csv("C:/Users/giuse/OneDrive/Desktop/Preparazione tesi/Thesis Code/OT Algorithm definitive/parameter_tuning_summary.csv")


#Top models
#1) 40 bins - 0.0050 eps (0.1279021)
#2) 75 bins - 0.0050 eps (0.1319141)
#3) 50 bins - 0.0050 eps (0.2624915)


# Standardized residuals
resids_df <- data.frame(
  ITA = resid_standardized_ITA,
  GER     = resid_standardized_GER,
  FRA     = resid_standardized_FRA,
  SPA    = resid_standardized_SPA
)

# ============================================ COMPARISON REAL DATA vs SIMULATED DATA ============================================
# ================================================================================================================================

# REAL DATA: RED
# SIMULATED DATA: GREEN

#--------------------------------------------- MARGINALS CHECK ------------------------------------------------------------------


index_labels <- c("ITA" = "FTSEMIB", "GER" = "DAX", "FRA" = "CAC", "SPA" = "IBEX")
vars <- c("ITA", "GER", "FRA", "SPA")

# ECDF COMPARISON
generate_ecdf_plots <- function(sim_data, real_data, n_sample = 1500) {
  plots <- list()
  set.seed(123)
  sim_sampled <- sim_data[sample(nrow(sim_data), n_sample), ]
  
  for (v in vars) {
    df_plot <- rbind(
      data.frame(value = real_data[[v]], type = "Real"),
      data.frame(value = sim_sampled[[v]], type = "Simulated")
    )
    
    p <- ggplot(df_plot, aes(x = value, color = type)) +
      stat_ecdf(geom = "step", size = 1) +
      scale_color_manual(values = c("Real" = "#e41a1c", "Simulated" = "#1b9e77")) +
      labs(title = index_labels[v], x = NULL, y = NULL) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 30, face = "bold", hjust = 0.5),  # Titoli più grandi
        axis.text = element_text(size = 9)
      )
    
    plots[[v]] <- p
  }
  return(plots)
}

for (i in 1:3) {
  sim_data <- sim_list[[i]]
  plots <- generate_ecdf_plots(sim_data, resids_df)
  
  grid.arrange(
    grobs = plots,
    ncol = 2
  )
}

for (i in 1:3) {
  sim_data <- sim_list[[i]]
  plots <- generate_ecdf_plots(sim_data, resids_df)
  
  ggsave(
    filename = paste0("Data/Chart Post preliminar analysis/ecdf_comparison_model_", i, ".pdf"),
    plot = arrangeGrob(
      grobs = plots,
      ncol = 2
    ),
    width = 10, height = 8
  )
}


# Ks test for ECDF
sim_sample_TOP1 <- sinkhorn_simulations_TOP_1[sample(nrow(sinkhorn_simulations_TOP_1), 1500), ]
sim_sample_TOP2 <- sinkhorn_simulations_TOP_2[sample(nrow(sinkhorn_simulations_TOP_2), 1500), ]
sim_sample_TOP3 <- sinkhorn_simulations_TOP_3[sample(nrow(sinkhorn_simulations_TOP_3), 1500), ]

ks.test(x=sim_sample_TOP1$ITA, y= resids_df$ITA, alternative = c("two.sided"))
ks.test(x=sim_sample_TOP1$GER, y= resids_df$GER, alternative = c("two.sided"))
ks.test(x=sim_sample_TOP1$FRA, y= resids_df$FRA, alternative = c("two.sided"))
ks.test(x=sim_sample_TOP1$SPA, y= resids_df$SPA, alternative = c("two.sided"))

ks.test(x=sim_sample_TOP2$ITA, y= resids_df$ITA, alternative = c("two.sided"))
ks.test(x=sim_sample_TOP2$GER, y= resids_df$GER, alternative = c("two.sided"))
ks.test(x=sim_sample_TOP2$FRA, y= resids_df$FRA, alternative = c("two.sided"))
ks.test(x=sim_sample_TOP2$SPA, y= resids_df$SPA, alternative = c("two.sided"))

ks.test(x=sim_sample_TOP3$ITA, y= resids_df$ITA, alternative = c("two.sided"))
ks.test(x=sim_sample_TOP3$GER, y= resids_df$GER, alternative = c("two.sided"))
ks.test(x=sim_sample_TOP3$FRA, y= resids_df$FRA, alternative = c("two.sided"))
ks.test(x=sim_sample_TOP3$SPA, y= resids_df$SPA, alternative = c("two.sided"))


#------------------------------------------ DEPENDENCE STRUCTURE CHECK ------------------------------------------------------------------

# Compare dependence sctructure (clouds of points)
plots_sinkhorn_dependence_resids <- function(sim_df, title_suffix = "") {
  combinations <- combn(c("ITA", "GER", "FRA", "SPA"), 2, simplify = FALSE)
  plots <- list()
  
  sim_sample <- sim_df[sample(nrow(sim_df), 1500), ]
  
  for (i in seq_along(combinations)) {
    pair <- combinations[[i]]
    var1 <- pair[1]
    var2 <- pair[2]
    
    x_real <- resids_df[[var1]]
    y_real <- resids_df[[var2]]
    x_sim  <- sim_sample[[var1]]
    y_sim  <- sim_sample[[var2]]
    
    df_real <- data.frame(x = x_real, y = y_real, type = "Real")
    df_sim  <- data.frame(x = x_sim,  y = y_sim,  type = "Simulated")
    df_plot <- rbind(df_real, df_sim)
    
    p <- ggplot(df_plot, aes(x = x, y = y, color = type)) +
      geom_point(alpha = 0.4, size = 1) +
      scale_color_manual(values = c("Real" = "red", "Simulated" = "#1b9e77")) +
      labs(title = NULL,
           x = index_labels[var1], y = index_labels[var2]) +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 8)
      )
    
    plots[[i]] <- p
  }
  
  return(plots)
}


plotsTOP1 <- plots_sinkhorn_dependence_resids(sinkhorn_simulations_TOP_1)
plotsTOP2 <- plots_sinkhorn_dependence_resids(sinkhorn_simulations_TOP_2)
plotsTOP3 <- plots_sinkhorn_dependence_resids(sinkhorn_simulations_TOP_3)

grob1 <- arrangeGrob(grobs = plotsTOP1, ncol = 3)
grob2 <- arrangeGrob(grobs = plotsTOP2, ncol = 3)
grob3 <- arrangeGrob(grobs = plotsTOP3, ncol = 3)

ggsave("Data/Chart Post preliminar analysis/dependence_sinkhorn_TOP1.pdf", grob1, width = 12, height = 8)
ggsave("Data/Chart Post preliminar analysis/dependence_sinkhorn_TOP2.pdf", grob2, width = 12, height = 8)
ggsave("Data/Chart Post preliminar analysis/dependence_sinkhorn_TOP3.pdf", grob3, width = 12, height = 8)

# ================================================================================================================================
# ================================================================================================================================


# ============================================== Algorithm for VaR ===============================================================
# ================================================================================================================================


# Parameters
window_size <- 1500  # Rolling window length
omega <- c(.25, .25, .25, .25)  # Portfolio weights

# Rolling window creation
rolling_ITA <- rollapply(X_tot$ITA, width = window_size, by = 1, FUN = identity, align = "right")
rolling_GER <- rollapply(X_tot$GER, width = window_size, by = 1, FUN = identity, align = "right")
rolling_FRA <- rollapply(X_tot$FRA, width = window_size, by = 1, FUN = identity, align = "right")
rolling_SPA  <- rollapply(X_tot$SPA,  width = window_size, by = 1, FUN = identity, align = "right")

num_forecasts <- nrow(rolling_ITA)  # Number of forecast points
test_indices <- 1:num_forecasts    # Index vector for parallel loop

# Parallel setup
num_cores <- parallel::detectCores() - 1
cl <- makeCluster(num_cores)

#  Set the stream for randmo numbers (RNG) for any worker
clusterSetRNGStream(cl, iseed = 123)


# Export required variables and functions
clusterExport(cl, varlist = c("rolling_ITA", "rolling_GER", "rolling_FRA", "rolling_SPA",
                              "spec_ITA", "spec_GER", "spec_FRA", "spec_SPA",
                              "simulate_logreturns", "uniform_trasformer", "omega",
                              "sinkhorn_simulations_TOP_1", "sinkhorn_simulations_TOP_2",
                              "sinkhorn_simulations_TOP_3"))

# Load packages in each worker
clusterEvalQ(cl, {
  library(rugarch)
  library(copula)
})

start_time <- Sys.time()

results_OT <- parLapply(cl, test_indices, function(i) {
  tryCatch({
    message("Step 1: Index i = ", i)
    
    sample_TOP_1 <- sinkhorn_simulations_TOP_1[sample(nrow(sinkhorn_simulations_TOP_1), 10000), ]
    sample_TOP_2 <- sinkhorn_simulations_TOP_2[sample(nrow(sinkhorn_simulations_TOP_2), 10000), ]
    sample_TOP_3 <- sinkhorn_simulations_TOP_3[sample(nrow(sinkhorn_simulations_TOP_3), 10000), ]
    
    window_data_ITA <- rolling_ITA[i, ]
    window_data_GER <- rolling_GER[i, ]
    window_data_FRA <- rolling_FRA[i, ]
    window_data_SPA <- rolling_SPA[i, ]
    message("Step 3: Rolling window extracted")
    
    fit_garch_ITA <- ugarchfit(spec = spec_ITA, data = window_data_ITA, solver = "hybrid")
    fit_garch_GER <- ugarchfit(spec = spec_GER, data = window_data_GER, solver = "hybrid")
    fit_garch_FRA <- ugarchfit(spec = spec_FRA, data = window_data_FRA, solver = "hybrid")
    fit_garch_SPA <- ugarchfit(spec = spec_SPA, data = window_data_SPA, solver = "hybrid")
    message("Step 4: GARCH fitted")
    
    if (any(c(fit_garch_ITA@fit$convergence, fit_garch_GER@fit$convergence,
              fit_garch_FRA@fit$convergence, fit_garch_SPA@fit$convergence) != 0)) {
      stop("At least one GARCH model did not converge")
    }
    message("Step 5: GARCH convergence OK")
    
    coef_ITA <- coef(fit_garch_ITA)
    coef_GER <- coef(fit_garch_GER)
    coef_FRA <- coef(fit_garch_FRA)
    coef_SPA <- coef(fit_garch_SPA)
    
    if (!all(c("shape", "skew") %in% names(coef_ITA))) stop("ITA shape/skew missing")
    if (!all(c("shape", "skew") %in% names(coef_GER))) stop("GER shape/skew missing")
    if (!all(c("shape", "skew") %in% names(coef_FRA))) stop("FRA shape/skew missing")
    if (!all(c("shape", "skew") %in% names(coef_SPA))) stop("SPA shape/skew missing")
    message("Step 6: Parameters OK")
    
    vol_forecast_ITA <- as.numeric(sigma(ugarchforecast(fit_garch_ITA, n.ahead = 1)))
    mean_forecast_ITA <- as.numeric(fitted(ugarchforecast(fit_garch_ITA, n.ahead = 1)))
    
    vol_forecast_GER <- as.numeric(sigma(ugarchforecast(fit_garch_GER, n.ahead = 1)))
    mean_forecast_GER <- as.numeric(fitted(ugarchforecast(fit_garch_GER, n.ahead = 1)))
    
    vol_forecast_FRA <- as.numeric(sigma(ugarchforecast(fit_garch_FRA, n.ahead = 1)))
    mean_forecast_FRA <- as.numeric(fitted(ugarchforecast(fit_garch_FRA, n.ahead = 1)))
    
    vol_forecast_SPA <- as.numeric(sigma(ugarchforecast(fit_garch_SPA, n.ahead = 1)))
    mean_forecast_SPA <- as.numeric(fitted(ugarchforecast(fit_garch_SPA, n.ahead = 1)))
    message("Step 7: Forecasts OK")
    
    logrets_TOP1_ITA <- mean_forecast_ITA + vol_forecast_ITA * sample_TOP_1$ITA
    logrets_TOP1_GER <- mean_forecast_GER + vol_forecast_GER *  sample_TOP_1$GER
    logrets_TOP1_FRA <- mean_forecast_FRA + vol_forecast_FRA *  sample_TOP_1$FRA
    logrets_TOP1_SPA <- mean_forecast_SPA + vol_forecast_SPA *  sample_TOP_1$SPA
    
    logrets_TOP2_ITA <- mean_forecast_ITA + vol_forecast_ITA * sample_TOP_2$ITA
    logrets_TOP2_GER <- mean_forecast_GER + vol_forecast_GER *  sample_TOP_2$GER
    logrets_TOP2_FRA <- mean_forecast_FRA + vol_forecast_FRA *  sample_TOP_2$FRA
    logrets_TOP2_SPA <- mean_forecast_SPA + vol_forecast_SPA *  sample_TOP_2$SPA
    
    logrets_TOP3_ITA <- mean_forecast_ITA + vol_forecast_ITA * sample_TOP_3$ITA
    logrets_TOP3_GER <- mean_forecast_GER + vol_forecast_GER *  sample_TOP_3$GER
    logrets_TOP3_FRA <- mean_forecast_FRA + vol_forecast_FRA *  sample_TOP_3$FRA
    logrets_TOP3_SPA <- mean_forecast_SPA + vol_forecast_SPA *  sample_TOP_3$SPA
    
    sim_logrets_TOP1 <- cbind(logrets_TOP1_ITA, logrets_TOP1_GER, logrets_TOP1_FRA, logrets_TOP1_SPA)
    sim_logrets_TOP2 <- cbind(logrets_TOP2_ITA, logrets_TOP2_GER, logrets_TOP2_FRA, logrets_TOP2_SPA)
    sim_logrets_TOP3 <- cbind(logrets_TOP3_ITA, logrets_TOP3_GER, logrets_TOP3_FRA, logrets_TOP3_SPA)
    
    sim_rets_TOP1 <- exp(sim_logrets_TOP1) - 1
    sim_rets_TOP2 <- exp(sim_logrets_TOP2) - 1
    sim_rets_TOP3 <- exp(sim_logrets_TOP3) - 1

    sim_rets_pl_TOP1 <- sim_rets_TOP1 %*% omega
    sim_rets_pl_TOP2 <- sim_rets_TOP2 %*% omega
    sim_rets_pl_TOP3 <- sim_rets_TOP3 %*% omega
    message("Step 9: Portfolio returns OK")
    
    list(
      Var95_OT_TOP1 = quantile(sim_rets_pl_TOP1, .05),
      VaR99_OT_TOP1 = quantile(sim_rets_pl_TOP1, .01),
      Var95_OT_TOP2 = quantile(sim_rets_pl_TOP2, .05),
      VaR99_OT_TOP2 = quantile(sim_rets_pl_TOP2, .01),
      Var95_OT_TOP3 = quantile(sim_rets_pl_TOP3, .05),
      VaR99_OT_TOP3 = quantile(sim_rets_pl_TOP3, .01)
    )
  }, error = function(e) {
    cat("Errore a i =", i, ":", conditionMessage(e), "\n")
    list(Var95_OT_TOP1 = NA, VaR99_OT_TOP1 = NA, Var95_OT_TOP2 = NA,
         VaR99_OT_TOP2 = NA, Var95_OT_TOP3 = NA, VaR99_OT_TOP3 = NA)
  })
})

stopCluster(cl)
end_time <- Sys.time()
cat("Time:", end_time - start_time, "\n")

# VaR extraction
Var95_OT_TOP1 <- sapply(results_OT, function(x) x$Var95_OT_TOP1)
VaR99_OT_TOP1 <- sapply(results_OT, function(x) x$VaR99_OT_TOP1)

Var95_OT_TOP2 <- sapply(results_OT, function(x) x$Var95_OT_TOP2)
VaR99_OT_TOP2 <- sapply(results_OT, function(x) x$VaR99_OT_TOP2)

Var95_OT_TOP3 <- sapply(results_OT, function(x) x$Var95_OT_TOP3)
VaR99_OT_TOP3 <- sapply(results_OT, function(x) x$VaR99_OT_TOP3)


Var95_OT_TOP1 <- Var95_OT_TOP1[1:1110]
VaR99_OT_TOP1 <- VaR99_OT_TOP1[1:1110]
Var95_OT_TOP2 <- Var95_OT_TOP2[1:1110]
VaR99_OT_TOP2 <- VaR99_OT_TOP2[1:1110]
Var95_OT_TOP3 <- Var95_OT_TOP3[1:1110]
VaR99_OT_TOP3 <- VaR99_OT_TOP3[1:1110]


# Save results
results_dir4 <- file.path("data", "results OPTIMAL TRANSPORT")
if (!dir.exists(results_dir4)) {
  dir.create(results_dir4, recursive = TRUE)
}
saveRDS(Var95_OT_TOP1, file.path(results_dir4, "Var95_OT_TOP1.rds"))
saveRDS(VaR99_OT_TOP1, file.path(results_dir4, "VaR99_OT_TOP1.rds"))
saveRDS(Var95_OT_TOP2, file.path(results_dir4, "Var95_OT_TOP2.rds"))
saveRDS(VaR99_OT_TOP2, file.path(results_dir4, "VaR99_OT_TOP2.rds"))
saveRDS(Var95_OT_TOP3, file.path(results_dir4, "Var95_OT_TOP3.rds"))
saveRDS(VaR99_OT_TOP3, file.path(results_dir4, "VaR99_OT_TOP3.rds"))


# Optimal Trasport plots
VaR95_list <- list(Var95_OT_TOP1, Var95_OT_TOP2, Var95_OT_TOP3)
VaR99_list <- list(VaR99_OT_TOP1, VaR99_OT_TOP2, VaR99_OT_TOP3)
titles <- c("40 bins - 0.0050 eps", "75 bins - 0.0050 eps", "50 bins - 0.0050 eps")

col_custom <- c("black", "darkred", "lightcoral", "red", "orange")
plots_list <- list()

for (i in 1:3) {
  p <- (function() {
    plot_VaR_OT(V_test = V_test[1:1110],
                dates = X_test$date,
                VaR95 = VaR95_list[[i]],
                VaR99 = VaR99_list[[i]],
                col = col_custom,
                title_text = titles[i])  
    ggplot2::last_plot()
  })()
  
  plots_list[[i]] <- p
}

combined_OT <- wrap_plots(plots_list, ncol = 1)

ggsave("Data/Chart Post preliminar analysis/combined_VaR_OT.pdf",
       combined_OT, width = 10, height = 12)


# Results

# Tests
test95_OT_TOP1 <- VaRTest(alpha = 0.05, actual = V_test, VaR = Var95_OT_TOP1)  
test99_OT_TOP1 <- VaRTest(alpha = 0.01, actual = V_test, VaR = VaR99_OT_TOP1)

test95_OT_TOP2 <- VaRTest(alpha = 0.05, actual = V_test, VaR = Var95_OT_TOP2)  
test99_OT_TOP2 <- VaRTest(alpha = 0.01, actual = V_test, VaR = VaR99_OT_TOP2)

test95_OT_TOP3 <- VaRTest(alpha = 0.05, actual = V_test, VaR = Var95_OT_TOP3)  
test99_OT_TOP3 <- VaRTest(alpha = 0.01, actual = V_test, VaR = VaR99_OT_TOP3)

ZT95_OT_TOP1 <- test95_OT_TOP1$actual.exceed / length(V_test)
ZT99_OT_TOP1 <- test99_OT_TOP1$actual.exceed / length(V_test)

ZT95_OT_TOP2 <- test95_OT_TOP2$actual.exceed / length(V_test)
ZT99_OT_TOP2 <- test99_OT_TOP2$actual.exceed / length(V_test)

ZT95_OT_TOP3 <- test95_OT_TOP3$actual.exceed / length(V_test)
ZT99_OT_TOP3 <- test99_OT_TOP3$actual.exceed / length(V_test)









models_OT <- c("40 bins-0.0050 eps", "75 bins-0.0050 eps", "50 bins-0.0050 eps")

# === 95% ===
ExpVio95_OT <- c(test95_OT_TOP1$expected.exceed,
                 test95_OT_TOP2$expected.exceed,
                 test95_OT_TOP3$expected.exceed)

ActVio95_OT <- c(test95_OT_TOP1$actual.exceed,
                 test95_OT_TOP2$actual.exceed,
                 test95_OT_TOP3$actual.exceed)

ZT95_OT <- round(ActVio95_OT / length(V_test), 4)

UCstat95_OT <- c(test95_OT_TOP1$uc.LRstat,
                 test95_OT_TOP2$uc.LRstat,
                 test95_OT_TOP3$uc.LRstat)

UCpval95_OT <- c(test95_OT_TOP1$uc.LRp,
                 test95_OT_TOP2$uc.LRp,
                 test95_OT_TOP3$uc.LRp)

CCstat95_OT <- c(test95_OT_TOP1$cc.LRstat,
                 test95_OT_TOP2$cc.LRstat,
                 test95_OT_TOP3$cc.LRstat)

CCpval95_OT <- c(test95_OT_TOP1$cc.LRp,
                 test95_OT_TOP2$cc.LRp,
                 test95_OT_TOP3$cc.LRp)

LRuc95_fmt <- paste0(round(UCstat95_OT, 2), " (", round(UCpval95_OT, 4), ")")
LRcc95_fmt <- paste0(round(CCstat95_OT, 2), " (", round(CCpval95_OT, 4), ")")

Lopez95_OT <- round(c(
  lopez_loss(V_test, Var95_OT_TOP1),
  lopez_loss(V_test, Var95_OT_TOP2),
  lopez_loss(V_test, Var95_OT_TOP3)
), 5)

# === Table  95% ===
tableVaR95_OT <- data.frame(
  Model = models_OT,
  ExpVio95 = ExpVio95_OT,
  ActVio95 = ActVio95_OT,
  MT95 = ZT95_OT,
  LRuc = LRuc95_fmt,
  LRcc = LRcc95_fmt,
  Lopez95 = Lopez95_OT * (100/ActVio95_OT)
)


# === 99% ===
ExpVio99_OT <- c(test99_OT_TOP1$expected.exceed,
                 test99_OT_TOP2$expected.exceed,
                 test99_OT_TOP3$expected.exceed)

ActVio99_OT <- c(test99_OT_TOP1$actual.exceed,
                 test99_OT_TOP2$actual.exceed,
                 test99_OT_TOP3$actual.exceed)

ZT99_OT <- round(ActVio99_OT / length(V_test), 4)

UCstat99_OT <- c(test99_OT_TOP1$uc.LRstat,
                 test99_OT_TOP2$uc.LRstat,
                 test99_OT_TOP3$uc.LRstat)

UCpval99_OT <- c(test99_OT_TOP1$uc.LRp,
                 test99_OT_TOP2$uc.LRp,
                 test99_OT_TOP3$uc.LRp)

CCstat99_OT <- c(test99_OT_TOP1$cc.LRstat,
                 test99_OT_TOP2$cc.LRstat,
                 test99_OT_TOP3$cc.LRstat)

CCpval99_OT <- c(test99_OT_TOP1$cc.LRp,
                 test99_OT_TOP2$cc.LRp,
                 test99_OT_TOP3$cc.LRp)

LRuc99_fmt <- paste0(round(UCstat99_OT, 2), " (", round(UCpval99_OT, 4), ")")
LRcc99_fmt <- paste0(round(CCstat99_OT, 2), " (", round(CCpval99_OT, 4), ")")

Lopez99_OT <- round(c(
  lopez_loss(V_test, VaR99_OT_TOP1),
  lopez_loss(V_test, VaR99_OT_TOP2),
  lopez_loss(V_test, VaR99_OT_TOP3)
), 5)

# === Table 99% ===
tableVaR99_OT <- data.frame(
  Model = models_OT,
  ExpVio99 = ExpVio99_OT,
  ActVio99 = ActVio99_OT,
  MT99 = ZT99_OT,
  LRuc = LRuc99_fmt,
  LRcc = LRcc99_fmt,
  Lopez99 = Lopez99_OT * (100/ActVio99_OT)
)