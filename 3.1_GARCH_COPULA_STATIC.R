# ==============================================================================
#                        GARCH-Copula Model with Static Fit
# ==============================================================================
# This script implements the algorithm used to compute the one-step-ahead 
# Value at Risk (VaR at T+1) based on a GARCH-Copula model with static copula 
# refitting (Copulas fitted only on training set). The procedure operates as follows:
#
# 1. Define rolling windows of 1500 observations for each time series using 
#    the `rollapply` function.
# 2. Fit the 6 Copulas on the training set

# 3. For each window, apply the best marginal model for each stock index, 
#    specifically ARMA(0,0) - TGARCH(1,1) with skewed t-distributed innovations.
# 4. Forecast the one-step-ahead conditional volatility for each index.
# 5. Transform the standardized residuals into standard uniform variables 
#    (pseudo-observations) via the Probability Integral Transform (PIT).
# 6. For each fitted copula, simulate 10,000 pseudo-observations via Monte Carlo.
# 7. Transform the simulated pseudo-observations into simulated log-returns by 
#    multiplying them by the forecasted volatilities.
# 8. Aggregate the simulated log-returns to obtain portfolio log-returns and 
#    convert them into arithmetic returns.
# 9. Estimate the Value at Risk by computing the 5th and 1st percentiles of the 
#    simulated distribution of portfolio returns.
#
# Steps 3–9 are repeated for each rolling window over the test set, which consists 
# of 1,110 observations.
#
# The entire algorithm is parallelized to improve computational efficiency.
# ==============================================================================



# STATIC COPULA FITTING (fitted only on the training set)
copula_static_t       <- fit.tCopula@copula
copula_static_gauss   <- fit.gaussCopula@copula
copula_static_frank   <- fit.frankCopula@copula
copula_static_gumbel  <- fit.gumbelCopula@copula
copula_static_clayton <- fit.claytonCopula@copula
copula_static_joe     <- fit.joeCopula@copula

# Parallel setup
num_cores <- parallel::detectCores() - 1
cl <- makeCluster(num_cores)

#  Set the stream for randmo numbers (RNG) for any worker
clusterSetRNGStream(cl, iseed = 123)

# Export required variables and functions
clusterExport(cl, varlist = c("rolling_ITA", "rolling_GER", "rolling_FRA", "rolling_SPA",
                              "spec_ITA", "spec_GER", "spec_FRA", "spec_SPA",
                              "simulate_logreturns", "uniform_trasformer", "omega",
                              "copula_static_t", "copula_static_gauss", "copula_static_frank",
                              "copula_static_gumbel", "copula_static_clayton", "copula_static_joe"))

# Load packages in each worker
clusterEvalQ(cl, {
  library(rugarch)
  library(copula)
})

start_time <- Sys.time()

# Rolling simulation with static copulas
results_static <- parLapply(cl, test_indices, function(i) {
  tryCatch({
    window_data_ITA <- rolling_ITA[i, ]
    window_data_GER <- rolling_GER[i, ]
    window_data_FRA <- rolling_FRA[i, ]
    window_data_SPA  <- rolling_SPA[i, ]
    
    fit_garch_ITA <- ugarchfit(spec = spec_ITA, data = window_data_ITA, solver = "hybrid")
    fit_garch_GER <- ugarchfit(spec = spec_GER, data = window_data_GER, solver = "hybrid")
    fit_garch_FRA <- ugarchfit(spec = spec_FRA, data = window_data_FRA, solver = "hybrid")
    fit_garch_SPA  <- ugarchfit(spec = spec_SPA,  data = window_data_SPA,  solver = "hybrid")
    
    rU_t       <- rCopula(10000, copula = copula_static_t)
    rU_g       <- rCopula(10000, copula = copula_static_gauss)
    rU_frank   <- rCopula(10000, copula = copula_static_frank)
    rU_gumbel  <- rCopula(10000, copula = copula_static_gumbel)
    rU_clayton <- rCopula(10000, copula = copula_static_clayton)
    rU_joe     <- rCopula(10000, copula = copula_static_joe)
    
    simulate_set <- function(fits, U) {
      logrets <- mapply(simulate_logreturns, fits, as.data.frame(U), SIMPLIFY = FALSE)
      cbind(exp(do.call(cbind, logrets)) - 1)
    }
    
    sim_rets_t  <- simulate_set(list(fit_garch_ITA, fit_garch_GER, fit_garch_FRA, fit_garch_SPA), rU_t)
    sim_rets_g  <- simulate_set(list(fit_garch_ITA, fit_garch_GER, fit_garch_FRA, fit_garch_SPA), rU_g)
    sim_rets_f  <- simulate_set(list(fit_garch_ITA, fit_garch_GER, fit_garch_FRA, fit_garch_SPA), rU_frank)
    sim_rets_gu <- simulate_set(list(fit_garch_ITA, fit_garch_GER, fit_garch_FRA, fit_garch_SPA), rU_gumbel)
    sim_rets_c  <- simulate_set(list(fit_garch_ITA, fit_garch_GER, fit_garch_FRA, fit_garch_SPA), rU_clayton)
    sim_rets_joe<- simulate_set(list(fit_garch_ITA, fit_garch_GER, fit_garch_FRA, fit_garch_SPA), rU_joe)
    
    Sim_PL <- function(sim_rets) sim_rets %*% omega
    
    list(
      VaR99_T_static = quantile(Sim_PL(sim_rets_t), 0.01),
      VaR95_T_static = quantile(Sim_PL(sim_rets_t), 0.05),
      VaR99_GAUSS_static = quantile(Sim_PL(sim_rets_g), 0.01),
      VaR95_GAUSS_static = quantile(Sim_PL(sim_rets_g), 0.05),
      VaR99_FRANK_static = quantile(Sim_PL(sim_rets_f), 0.01),
      VaR95_FRANK_static = quantile(Sim_PL(sim_rets_f), 0.05),
      VaR99_GUMBEL_static = quantile(Sim_PL(sim_rets_gu), 0.01),
      VaR95_GUMBEL_static = quantile(Sim_PL(sim_rets_gu), 0.05),
      VaR99_CLAYTON_static = quantile(Sim_PL(sim_rets_c), 0.01),
      VaR95_CLAYTON_static = quantile(Sim_PL(sim_rets_c), 0.05),
      VaR99_JOE_static = quantile(Sim_PL(sim_rets_joe), 0.01),
      VaR95_JOE_static = quantile(Sim_PL(sim_rets_joe), 0.05)
    )
  }, error = function(e) {
    list(
      VaR99_T_static = NA, VaR95_T_static = NA,
      VaR99_GAUSS_static = NA, VaR95_GAUSS_static = NA,
      VaR99_FRANK_static = NA, VaR95_FRANK_static = NA,
      VaR99_GUMBEL_static = NA, VaR95_GUMBEL_static = NA,
      VaR99_CLAYTON_static = NA, VaR95_CLAYTON_static = NA,
      VaR99_JOE_static = NA, VaR95_JOE_static = NA
    )
  })
})

stopCluster(cl)
end_time <- Sys.time()
cat("Time:", end_time - start_time, "\n")

# Extract static copula VaRs
VaR99_MC_T_static <- sapply(results_static, function(x) x$VaR99_T_static)
VaR95_MC_T_static <- sapply(results_static, function(x) x$VaR95_T_static)
VaR99_MC_GAUSS_static <- sapply(results_static, function(x) x$VaR99_GAUSS_static)
VaR95_MC_GAUSS_static <- sapply(results_static, function(x) x$VaR95_GAUSS_static)
VaR99_MC_FRANK_static <- sapply(results_static, function(x) x$VaR99_FRANK_static)
VaR95_MC_FRANK_static <- sapply(results_static, function(x) x$VaR95_FRANK_static)
VaR99_MC_GUMBEL_static <- sapply(results_static, function(x) x$VaR99_GUMBEL_static)
VaR95_MC_GUMBEL_static <- sapply(results_static, function(x) x$VaR95_GUMBEL_static)
VaR99_MC_CLAYTON_static <- sapply(results_static, function(x) x$VaR99_CLAYTON_static)
VaR95_MC_CLAYTON_static <- sapply(results_static, function(x) x$VaR95_CLAYTON_static)
VaR99_MC_JOE_static <- sapply(results_static, function(x) x$VaR99_JOE_static)
VaR95_MC_JOE_static <- sapply(results_static, function(x) x$VaR95_JOE_static)

VaR99_MC_T_static <- VaR99_MC_T_static[1:1110]
VaR95_MC_T_static <- VaR95_MC_T_static[1:1110]
VaR99_MC_GAUSS_static <- VaR99_MC_GAUSS_static[1:1110]
VaR95_MC_GAUSS_static <- VaR95_MC_GAUSS_static[1:1110]
VaR99_MC_FRANK_static <- VaR99_MC_FRANK_static[1:1110]
VaR95_MC_FRANK_static <- VaR95_MC_FRANK_static[1:1110]
VaR99_MC_GUMBEL_static <- VaR99_MC_GUMBEL_static[1:1110]
VaR95_MC_GUMBEL_static <- VaR95_MC_GUMBEL_static[1:1110]
VaR99_MC_CLAYTON_static <- VaR99_MC_CLAYTON_static[1:1110]
VaR95_MC_CLAYTON_static <- VaR95_MC_CLAYTON_static[1:1110]
VaR99_MC_JOE_static <- VaR99_MC_JOE_static[1:1110]
VaR95_MC_JOE_static <- VaR95_MC_JOE_static[1:1110]

# Results saving

results_dir2 <- file.path("data", "GARCH_COPULA STATIC")
if (!dir.exists(results_dir2)) {
  dir.create(results_dir2, recursive = TRUE)
}

saveRDS(VaR99_MC_T_static,        file.path(results_dir2, "VaR99_MC_tcopula_static.rds"))
saveRDS(VaR95_MC_T_static,        file.path(results_dir2, "VaR95_MC_tcopula_static.rds"))
saveRDS(VaR99_MC_GAUSS_static,    file.path(results_dir2, "VaR99_MC_gausscopula_static.rds"))
saveRDS(VaR95_MC_GAUSS_static,    file.path(results_dir2, "VaR95_MC_gausscopula_static.rds"))
saveRDS(VaR99_MC_FRANK_static,    file.path(results_dir2, "VaR99_MC_frankcopula_static.rds"))
saveRDS(VaR95_MC_FRANK_static,    file.path(results_dir2, "VaR95_MC_frankcopula_static.rds"))
saveRDS(VaR99_MC_GUMBEL_static,   file.path(results_dir2, "VaR99_MC_gumbelcopula_static.rds"))
saveRDS(VaR95_MC_GUMBEL_static,   file.path(results_dir2, "VaR95_MC_gumbelcopula_static.rds"))
saveRDS(VaR99_MC_CLAYTON_static,  file.path(results_dir2, "VaR99_MC_claytoncopula_static.rds"))
saveRDS(VaR95_MC_CLAYTON_static,  file.path(results_dir2, "VaR95_MC_claytoncopula_static.rds"))
saveRDS(VaR99_MC_JOE_static,      file.path(results_dir2, "VaR99_MC_joecopula_static.rds"))
saveRDS(VaR95_MC_JOE_static,      file.path(results_dir2, "VaR95_MC_joecopula_static.rds"))


#------------------------- Tests -------------------------------------------

test95_GAUSS_static <- VaRTest(alpha = 0.05, actual = V_test, VaR = VaR95_MC_GAUSS_static)  
test99_GAUSS_static <- VaRTest(alpha = 0.01, actual = V_test, VaR =VaR99_MC_GAUSS_static)

test95_T_static <- VaRTest(alpha = 0.05, actual = V_test, VaR = VaR95_MC_T_static)  
test99_T_static <- VaRTest(alpha = 0.01, actual = V_test, VaR =VaR99_MC_T_static)

test95_GUMBEL_static <- VaRTest(alpha = 0.05, actual = V_test, VaR = VaR95_MC_GUMBEL_static)  
test99_GUMBEL_static <- VaRTest(alpha = 0.01, actual = V_test, VaR =VaR99_MC_GUMBEL_static)

test95_CLAYTON_static <- VaRTest(alpha = 0.05, actual = V_test, VaR = VaR95_MC_CLAYTON_static)  
test99_CLAYTON_static <- VaRTest(alpha = 0.01, actual = V_test, VaR =VaR99_MC_CLAYTON_static)

test95_FRANK_static <- VaRTest(alpha = 0.05, actual = V_test, VaR = VaR95_MC_FRANK_static)  
test99_FRANK_static <- VaRTest(alpha = 0.01, actual = V_test, VaR =VaR99_MC_FRANK_static)

test95_JOE_static <- VaRTest(alpha = 0.05, actual = V_test, VaR = VaR95_MC_JOE_static)  
test99_JOE_static <- VaRTest(alpha = 0.01, actual = V_test, VaR =VaR99_MC_JOE_static)




# GARCH-Copula with statif fit charts
copula_names <- c("T", "GAUSS", "FRANK", "GUMBEL", "CLAYTON", "JOE")
plot_list_static <- vector("list", length(copula_names))

for (i in seq_along(copula_names)) {
  copula <- copula_names[i]
  VaR95 <- get(paste0("VaR95_MC_", copula, "_static"))
  VaR99 <- get(paste0("VaR99_MC_", copula, "_static"))
  
  p <- ggplot2::last_plot()  # placeholder
  p <- (function() {
    plot_VaR_hits(copula, V_test, X_test, VaR95, VaR99,
                  col = c("black", "darkgreen", "green", "red", "orange"))
    ggplot2::last_plot()
  })()
  
  p <- p + labs(title = copula, x = NULL, y = NULL) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 30),
          axis.text = element_text(size = 6))
  
  plot_list_static[[i]] <- p
}

combined_COPULA_STATIC <- wrap_plots(plot_list_static, ncol = 2, nrow = 3)

ggsave("Data/Chart Post preliminar analysis/STATIC_GARCH_COPULA.pdf", combined_COPULA_STATIC, width = 12, height = 8)



# ------------------------------------------- Table ---------------------------------------------


# Lopez loss function
# t student
lopez95_T_static <- lopez_loss(V_test, VaR95_MC_T_static)
lopez99_T_static <- lopez_loss(V_test, VaR99_MC_T_static)

# Gauss
lopez99_GAUSS_static<- lopez_loss(V_test, VaR99_MC_GAUSS_static)
lopez95_GAUSS_static<- lopez_loss(V_test, VaR95_MC_GAUSS_static)

# Clayton
lopez99_CLAYTON_static<- lopez_loss(V_test, VaR99_MC_CLAYTON_static)
lopez95_CLAYTON_static <- lopez_loss(V_test, VaR95_MC_CLAYTON_static)

# FRANK
lopez99_FRANK_static<- lopez_loss(V_test, VaR99_MC_FRANK_static)
lopez95_FRANK_static <- lopez_loss(V_test, VaR95_MC_FRANK_static)

# GUMBEL
lopez99_GUMBEL_static<- lopez_loss(V_test, VaR99_MC_GUMBEL_static)
lopez95_GUMBEL_static <- lopez_loss(V_test, VaR95_MC_GUMBEL_static)

# JOE
lopez99_JOE_static<- lopez_loss(V_test, VaR99_MC_JOE_static)
lopez95_JOE_static <- lopez_loss(V_test, VaR95_MC_JOE_static)

loss_df_static <- data.frame(
  Copula = c("t", "Gaussian", "Clayton", "Frank", "Gumbel", "Joe"),
  Lopez95 = c(lopez95_T_static, lopez95_GAUSS_static, lopez95_CLAYTON_static, lopez95_FRANK_static, lopez95_GUMBEL_static, lopez95_JOE_static),
  Lopez99 = c(lopez99_T_static,lopez99_GAUSS_static, lopez99_CLAYTON_static, lopez99_FRANK_static, lopez99_GUMBEL_static, lopez99_JOE_static)
  
)





copula_names_final <- c("Gaussian", "t", "Clayton", "Frank", "Gumbel", "Joe")

# === Combine stats e p-value ===
combine_stat_pval <- function(stat, pval) {
  paste0(round(stat, 2), " (", round(pval, 4), ")")
}

# ===  Table 95% ===
ExpVio95_static <- c(
  test95_GAUSS_static$expected.exceed,
  test95_T_static$expected.exceed,
  test95_CLAYTON_static$expected.exceed,
  test95_FRANK_static$expected.exceed,
  test95_GUMBEL_static$expected.exceed,
  test95_JOE_static$expected.exceed
)

ActVio95_static <- c(
  test95_GAUSS_static$actual.exceed,
  test95_T_static$actual.exceed,
  test95_CLAYTON_static$actual.exceed,
  test95_FRANK_static$actual.exceed,
  test95_GUMBEL_static$actual.exceed,
  test95_JOE_static$actual.exceed
)

ZT95_static <- ActVio95_static / length(V_test)

LRuc_95_stat <- c(
  test95_GAUSS_static$uc.LRstat,
  test95_T_static$uc.LRstat,
  test95_CLAYTON_static$uc.LRstat,
  test95_FRANK_static$uc.LRstat,
  test95_GUMBEL_static$uc.LRstat,
  test95_JOE_static$uc.LRstat
)

LRuc_95_pval_stat <- c(
  test95_GAUSS_static$uc.LRp,
  test95_T_static$uc.LRp,
  test95_CLAYTON_static$uc.LRp,
  test95_FRANK_static$uc.LRp,
  test95_GUMBEL_static$uc.LRp,
  test95_JOE_static$uc.LRp
)

LRcc_95_stat <- c(
  test95_GAUSS_static$cc.LRstat,
  test95_T_static$cc.LRstat,
  test95_CLAYTON_static$cc.LRstat,
  test95_FRANK_static$cc.LRstat,
  test95_GUMBEL_static$cc.LRstat,
  test95_JOE_static$cc.LRstat
)

LRcc_95_pval_stat <- c(
  test95_GAUSS_static$cc.LRp,
  test95_T_static$cc.LRp,
  test95_CLAYTON_static$cc.LRp,
  test95_FRANK_static$cc.LRp,
  test95_GUMBEL_static$cc.LRp,
  test95_JOE_static$cc.LRp
)

Lopez95_static <- c(
  lopez95_GAUSS_static,
  lopez95_T_static,
  lopez95_CLAYTON_static,
  lopez95_FRANK_static,
  lopez95_GUMBEL_static,
  lopez95_JOE_static
)

tableVaR95_static <- data.frame(
  Model = copula_names_final,
  ExpVio95 = ExpVio95_static,
  ActVio95 = ActVio95_static,
  MT95 = round(ZT95_static, 4),
  LRuc = combine_stat_pval(LRuc_95_stat, LRuc_95_pval_stat),
  LRcc = combine_stat_pval(LRcc_95_stat, LRcc_95_pval_stat),
  Lopez95 = round(Lopez95_static, 5) * (100/ActVio95_static)
)

# === Table 99% ===
ExpVio99_static <- c(
  test99_GAUSS_static$expected.exceed,
  test99_T_static$expected.exceed,
  test99_CLAYTON_static$expected.exceed,
  test99_FRANK_static$expected.exceed,
  test99_GUMBEL_static$expected.exceed,
  test99_JOE_static$expected.exceed
)

ActVio99_static <- c(
  test99_GAUSS_static$actual.exceed,
  test99_T_static$actual.exceed,
  test99_CLAYTON_static$actual.exceed,
  test99_FRANK_static$actual.exceed,
  test99_GUMBEL_static$actual.exceed,
  test99_JOE_static$actual.exceed
)

ZT99_static <- ActVio99_static / length(V_test)

LRuc_99_stat <- c(
  test99_GAUSS_static$uc.LRstat,
  test99_T_static$uc.LRstat,
  test99_CLAYTON_static$uc.LRstat,
  test99_FRANK_static$uc.LRstat,
  test99_GUMBEL_static$uc.LRstat,
  test99_JOE_static$uc.LRstat
)

LRuc_99_pval_stat <- c(
  test99_GAUSS_static$uc.LRp,
  test99_T_static$uc.LRp,
  test99_CLAYTON_static$uc.LRp,
  test99_FRANK_static$uc.LRp,
  test99_GUMBEL_static$uc.LRp,
  test99_JOE_static$uc.LRp
)

LRcc_99_stat <- c(
  test99_GAUSS_static$cc.LRstat,
  test99_T_static$cc.LRstat,
  test99_CLAYTON_static$cc.LRstat,
  test99_FRANK_static$cc.LRstat,
  test99_GUMBEL_static$cc.LRstat,
  test99_JOE_static$cc.LRstat
)

LRcc_99_pval_stat <- c(
  test99_GAUSS_static$cc.LRp,
  test99_T_static$cc.LRp,
  test99_CLAYTON_static$cc.LRp,
  test99_FRANK_static$cc.LRp,
  test99_GUMBEL_static$cc.LRp,
  test99_JOE_static$cc.LRp
)

Lopez99_static <- c(
  lopez99_GAUSS_static,
  lopez99_T_static,
  lopez99_CLAYTON_static,
  lopez99_FRANK_static,
  lopez99_GUMBEL_static,
  lopez99_JOE_static
)

# Table 99%
tableVaR99_static <- data.frame(
  Model = copula_names_final,
  ExpVio99 = ExpVio99_static,
  ActVio99 = ActVio99_static,
  MT99 = round(ZT99_static, 4),
  LRuc = combine_stat_pval(LRuc_99_stat, LRuc_99_pval_stat),
  LRcc = combine_stat_pval(LRcc_99_stat, LRcc_99_pval_stat),
  Lopez99 = round(Lopez99_static, 5) * (100/ActVio99_static)
)

