# ==============================================================================
#                        GARCH-Copula Model with Dynamic Fit
# ==============================================================================
# This script implements the algorithm used to compute the one-step-ahead 
# Value at Risk (VaR at T+1) based on a GARCH-Copula model with dynamic copula 
# refitting. The procedure operates as follows:
#
# 1. Define rolling windows of 1500 observations for each time series using 
#    the `rollapply` function.
# 2. For each window, apply the best marginal model for each stock index, 
#    specifically ARMA(0,0) - TGARCH(1,1) with skewed t-distributed innovations.
# 3. Forecast the one-step-ahead conditional volatility for each index.
# 4. Transform the standardized residuals into standard uniform variables 
#    (pseudo-observations) via the Probability Integral Transform (PIT).
# 5. Fit six different copula models to the pseudo-observations.
# 6. For each fitted copula, simulate 10,000 pseudo-observations via Monte Carlo.
# 7. Transform the simulated pseudo-observations into simulated log-returns by 
#    multiplying them by the forecasted volatilities.
# 8. Aggregate the simulated log-returns to obtain portfolio log-returns and 
#    convert them into arithmetic returns.
# 9. Estimate the Value at Risk by computing the 5th and 1st percentiles of the 
#    simulated distribution of portfolio returns.
#
# Steps 2–9 are repeated for each rolling window over the test set, which consists 
# of 1,110 observations.
#
# The entire algorithm is parallelized to improve computational efficiency.
#
#
# At the end of the script, a parameter stability analysis is performed to assess
# whether the marginal model parameters and copula parameters remain stable over time.
# ==============================================================================


# Parameters
window_size <- 1500  # Rolling window length
omega <- c(.25, .25, .25, .25)  # Portfolio weights

# Rolling window creation
rolling_ITA <- rollapply(X_tot$ITA, width = window_size, by = 1, FUN = identity, align = "right")
rolling_GER <- rollapply(X_tot$GER, width = window_size, by = 1, FUN = identity, align = "right")
rolling_FRA <- rollapply(X_tot$FRA, width = window_size, by = 1, FUN = identity, align = "right")
rolling_SPA  <- rollapply(X_tot$SPA,  width = window_size, by = 1, FUN = identity, align = "right")

num_forecasts <- nrow(rolling_ITA)  # Number of forecast points
test_indices <- 1:num_forecasts     # Index vector for parallel loop

# Parallel setup
num_cores <- parallel::detectCores() - 1  # Use all available cores minus one
cl <- makeCluster(num_cores)

#  Set the stream for randmo numbers (RNG) for any worker
clusterSetRNGStream(cl, iseed = 123)

# Export required variables and functions
clusterExport(cl, varlist = c("rolling_ITA", "rolling_GER", "rolling_FRA", "rolling_SPA",
                              "spec_ITA", "spec_GER", "spec_FRA", "spec_SPA",
                              "simulate_logreturns", "uniform_trasformer", "omega"))

# Load packages in each worker
clusterEvalQ(cl, {
  library(rugarch)
  library(copula)
})

start_time <- Sys.time()

# Rolling estimation and simulation
results_base <- parLapply(cl, test_indices, function(i) {
  tryCatch({
    # Extract data window
    window_data_ITA <- rolling_ITA[i, ]
    window_data_GER <- rolling_GER[i, ]
    window_data_FRA <- rolling_FRA[i, ]
    window_data_SPA  <- rolling_SPA[i, ]
    
    # Fit GARCH models for each series
    fit_garch_ITA <- ugarchfit(spec_ITA, data = window_data_ITA, solver = "hybrid")
    fit_garch_GER <- ugarchfit(spec_GER, data = window_data_GER, solver = "hybrid")
    fit_garch_FRA <- ugarchfit(spec_FRA, data = window_data_FRA, solver = "hybrid")
    fit_garch_SPA  <- ugarchfit(spec_SPA,  data = window_data_SPA,  solver = "hybrid")
    
    # Save GARCH parameters
    garch_params <- list(
      ITA = coef(fit_garch_ITA),
      GER = coef(fit_garch_GER),
      FRA = coef(fit_garch_FRA),
      SPA = coef(fit_garch_SPA)
    )
    
    # Transform residuals into PITs
    unif_matr <- uniform_trasformer(list(fit_garch_ITA, fit_garch_GER, fit_garch_FRA, fit_garch_SPA))
    
    # Fit copulas to PIT matrix
    fit_t_copula      <- fitCopula(tCopula(dim = 4, dispstr = "un"), data = unif_matr)
    fit_gauss_copula  <- fitCopula(normalCopula(dim = 4, dispstr = "un"), data = unif_matr)
    fit_frank_copula  <- fitCopula(frankCopula(dim = 4), data = unif_matr)
    fit_gumbel_copula <- fitCopula(gumbelCopula(dim = 4), data = unif_matr)
    fit_clayton_copula<- fitCopula(claytonCopula(dim = 4), data = unif_matr)
    fit_joe_copula    <- fitCopula(joeCopula(dim = 4), data = unif_matr)
    
    # Simulate PITs from each copula
    rU_t       <- rCopula(10000, copula = fit_t_copula@copula)
    rU_g       <- rCopula(10000, copula = fit_gauss_copula@copula)
    rU_frank   <- rCopula(10000, copula = fit_frank_copula@copula)
    rU_gumbel  <- rCopula(10000, copula = fit_gumbel_copula@copula)
    rU_clayton <- rCopula(10000, copula = fit_clayton_copula@copula)
    rU_joe     <- rCopula(10000, copula = fit_joe_copula@copula)
    
    # Function to simulate log-returns and convert to simple returns
    simulate_set <- function(fits, U) {
      logrets <- mapply(simulate_logreturns, fits, as.data.frame(U), SIMPLIFY = FALSE)
      cbind(exp(do.call(cbind, logrets)) - 1)
    }
    

    
    # Generate simulated returns
    sim_rets_t  <- simulate_set(list(fit_garch_ITA, fit_garch_GER, fit_garch_FRA, fit_garch_SPA), rU_t)
    sim_rets_g  <- simulate_set(list(fit_garch_ITA, fit_garch_GER, fit_garch_FRA, fit_garch_SPA), rU_g)
    sim_rets_f  <- simulate_set(list(fit_garch_ITA, fit_garch_GER, fit_garch_FRA, fit_garch_SPA), rU_frank)
    sim_rets_gu <- simulate_set(list(fit_garch_ITA, fit_garch_GER, fit_garch_FRA, fit_garch_SPA), rU_gumbel)
    sim_rets_c  <- simulate_set(list(fit_garch_ITA, fit_garch_GER, fit_garch_FRA, fit_garch_SPA), rU_clayton)
    sim_rets_joe<- simulate_set(list(fit_garch_ITA, fit_garch_GER, fit_garch_FRA, fit_garch_SPA), rU_joe)
    
    # Portfolio P&L from simulated returns
    Sim_PL <- function(sim_rets) sim_rets %*% omega
    
    # Return VaR estimates and copula parameters
    list(
      VaR99_T = quantile(Sim_PL(sim_rets_t), 0.01),
      VaR95_T = quantile(Sim_PL(sim_rets_t), 0.05),
      VaR99_GAUSS = quantile(Sim_PL(sim_rets_g), 0.01),
      VaR95_GAUSS = quantile(Sim_PL(sim_rets_g), 0.05),
      VaR99_FRANK = quantile(Sim_PL(sim_rets_f), 0.01),
      VaR95_FRANK = quantile(Sim_PL(sim_rets_f), 0.05),
      VaR99_GUMBEL = quantile(Sim_PL(sim_rets_gu), 0.01),
      VaR95_GUMBEL = quantile(Sim_PL(sim_rets_gu), 0.05),
      VaR99_CLAYTON = quantile(Sim_PL(sim_rets_c), 0.01),
      VaR95_CLAYTON = quantile(Sim_PL(sim_rets_c), 0.05),
      VaR99_JOE = quantile(Sim_PL(sim_rets_joe), 0.01),
      VaR95_JOE = quantile(Sim_PL(sim_rets_joe), 0.05),
      copula_params_t = coef(fit_t_copula),
      copula_params_g = coef(fit_gauss_copula),
      copula_params_frank = coef(fit_frank_copula),
      copula_params_gumbel = coef(fit_gumbel_copula),
      copula_params_clayton = coef(fit_clayton_copula),
      copula_params_joe = coef(fit_joe_copula),
      garch_params = garch_params
    )
  }, error = function(e) {
    # In case of error, return NA-filled result
    list(
      VaR99_T = NA, VaR95_T = NA, VaR99_GAUSS = NA, VaR95_GAUSS = NA,
      VaR99_FRANK = NA, VaR95_FRANK = NA, VaR99_GUMBEL = NA, VaR95_GUMBEL = NA,
      VaR99_CLAYTON = NA, VaR95_CLAYTON = NA, VaR99_JOE = NA, VaR95_JOE = NA,
      copula_params_t = rep(NA, 7), copula_params_g = rep(NA, 6),
      copula_params_frank = NA, copula_params_gumbel = NA,
      copula_params_clayton = NA, copula_params_joe = NA,
      garch_params = list(ITA = NA, GER = NA, FRA = NA, SPA = NA)
    )
  })
})

stopCluster(cl)  # Close cluster
cat("Time:", Sys.time() - start_time, "\n")  # Print execution time

# results extraction
VaR99_MC_T <- sapply(results_base, function(x) x$VaR99_T)
VaR95_MC_T <- sapply(results_base, function(x) x$VaR95_T)
VaR99_MC_GAUSS <- sapply(results_base, function(x) x$VaR99_GAUSS)
VaR95_MC_GAUSS <- sapply(results_base, function(x) x$VaR95_GAUSS)
VaR99_MC_FRANK <- sapply(results_base, function(x) x$VaR99_FRANK)
VaR95_MC_FRANK <- sapply(results_base, function(x) x$VaR95_FRANK)
VaR99_MC_GUMBEL <- sapply(results_base, function(x) x$VaR99_GUMBEL)
VaR95_MC_GUMBEL <- sapply(results_base, function(x) x$VaR95_GUMBEL)
VaR99_MC_CLAYTON <- sapply(results_base, function(x) x$VaR99_CLAYTON)
VaR95_MC_CLAYTON <- sapply(results_base, function(x) x$VaR95_CLAYTON)
VaR99_MC_JOE <- sapply(results_base, function(x) x$VaR99_JOE)
VaR95_MC_JOE <- sapply(results_base, function(x) x$VaR95_JOE)

# Excluding the last value because it has no match with the test portfolio
VaR99_MC_T <- VaR99_MC_T[1:1110]
VaR95_MC_T <- VaR95_MC_T[1:1110]
VaR99_MC_GAUSS <-  VaR99_MC_GAUSS[1:1110]
VaR95_MC_GAUSS <- VaR95_MC_GAUSS[1:1110]
VaR99_MC_FRANK <- VaR99_MC_FRANK[1:1110]
VaR95_MC_FRANK <- VaR95_MC_FRANK[1:1110]
VaR99_MC_GUMBEL <-VaR99_MC_GUMBEL[1:1110]
VaR95_MC_GUMBEL <- VaR95_MC_GUMBEL[1:1110]
VaR99_MC_CLAYTON <- VaR99_MC_CLAYTON[1:1110]
VaR95_MC_CLAYTON <- VaR95_MC_CLAYTON[1:1110]
VaR99_MC_JOE <- VaR99_MC_JOE[1:1110]
VaR95_MC_JOE <- VaR95_MC_JOE[1:1110]


copula_param_df_t <- as.data.frame(do.call(rbind, lapply(results_base, function(x) x$copula_params_t)))
copula_param_df_gauss <- as.data.frame(do.call(rbind, lapply(results_base, function(x) x$copula_params_g)))
copula_param_df_frank <- as.data.frame(do.call(rbind, lapply(results_base, function(x) x$copula_params_frank)))
copula_param_df_gumbel <- as.data.frame(do.call(rbind, lapply(results_base, function(x) x$copula_params_gumbel)))
copula_param_df_clayton <- as.data.frame(do.call(rbind, lapply(results_base, function(x) x$copula_params_clayton)))
copula_param_df_joe <- as.data.frame(do.call(rbind, lapply(results_base, function(x) x$copula_params_joe)))

copula_param_df_joe <- as.data.frame(do.call(rbind, lapply(results_base, function(x) x$copula_params_joe)))
garch_param_df <- lapply(results_base, function(x) x$garch_params)

nrow(copula_param_df_t)
nrow(X_test)

# Savings results
results_dir <- file.path("Data", "results GARCH_COPULA CLASSIC")
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

saveRDS(VaR99_MC_T,        file.path(results_dir, "VaR99_MC_tcopula.rds"))
saveRDS(VaR95_MC_T,        file.path(results_dir, "VaR95_MC_tcopula.rds"))
saveRDS(VaR99_MC_GAUSS,    file.path(results_dir, "VaR99_MC_gausscopula.rds"))
saveRDS(VaR95_MC_GAUSS,    file.path(results_dir, "VaR95_MC_gausscopula.rds"))
saveRDS(VaR99_MC_FRANK,    file.path(results_dir, "VaR99_MC_frankcopula.rds"))
saveRDS(VaR95_MC_FRANK,    file.path(results_dir, "VaR95_MC_frankcopula.rds"))
saveRDS(VaR99_MC_GUMBEL,   file.path(results_dir, "VaR99_MC_gumbelcopula.rds"))
saveRDS(VaR95_MC_GUMBEL,   file.path(results_dir, "VaR95_MC_gumbelcopula.rds"))
saveRDS(VaR99_MC_CLAYTON,  file.path(results_dir, "VaR99_MC_claytoncopula.rds"))
saveRDS(VaR95_MC_CLAYTON,  file.path(results_dir, "VaR95_MC_claytoncopula.rds"))
saveRDS(VaR99_MC_JOE,      file.path(results_dir, "VaR99_MC_joecopula.rds"))
saveRDS(VaR95_MC_JOE,      file.path(results_dir, "VaR95_MC_joecopula.rds"))

saveRDS(copula_param_df_t,       file.path(results_dir, "copula_params_tcopula.rds"))
saveRDS(copula_param_df_gauss,   file.path(results_dir, "copula_params_gausscopula.rds"))
saveRDS(copula_param_df_frank,   file.path(results_dir, "copula_params_frankcopula.rds"))
saveRDS(copula_param_df_gumbel,  file.path(results_dir, "copula_params_gumbelcopula.rds"))
saveRDS(copula_param_df_clayton, file.path(results_dir, "copula_params_claytoncopula.rds"))
saveRDS(copula_param_df_joe,     file.path(results_dir, "copula_params_joecopula.rds"))

saveRDS(garch_param_df, file.path(results_dir, "garch_parameters_timeseries.rds"))


V_test <- as.matrix(exp(as.matrix(X_test[,2:5])) - 1) %*% as.numeric(omega)


nrow(X_test["date"])
nrow(V_test)
length(VaR99_MC_T)
length(VaR99_MC_GAUSS)
length(VaR95_MC_T)
length(VaR95_MC_GAUSS)



# GARCH-Copula with dynamic fit CHARTS
copula_names <- c("T", "GAUSS", "FRANK", "GUMBEL", "CLAYTON", "JOE")
plot_list_dynamic <- vector("list", length(copula_names))

# Generate the six plots and store them
for (i in seq_along(copula_names)) {
  copula <- copula_names[i]
  VaR95 <- get(paste0("VaR95_MC_", copula))
  VaR99 <- get(paste0("VaR99_MC_", copula))
  
  p <- (function() {
    plot_VaR_hits(copula, V_test, X_test, VaR95, VaR99,
                  col = c("black", "darkblue", "lightblue", "red", "orange"))
    ggplot2::last_plot()
  })()
  
  p <- p + labs(title = copula, x = NULL, y = NULL) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          axis.text = element_text(size = 6))
  
  plot_list_dynamic[[i]] <- p  
}

combined_COPULA_DYNAMIC <- wrap_plots(plot_list_dynamic, ncol = 2, nrow = 3)

ggsave("Data/Chart Post preliminar analysis/DYNAMIC_GARCH_COPULA.pdf", 
       combined_COPULA_DYNAMIC, width = 12, height = 8)




#---------------------- Unconditional Test and Conditional Test ----------------------------------------
test95_GAUSS <- VaRTest(alpha = 0.05, actual = V_test, VaR = VaR95_MC_GAUSS)  
test99_GAUSS <- VaRTest(alpha = 0.01, actual = V_test, VaR =VaR99_MC_GAUSS)

test95_T <- VaRTest(alpha = 0.05, actual = V_test, VaR = VaR95_MC_T)  
test99_T <- VaRTest(alpha = 0.01, actual = V_test, VaR =VaR99_MC_T)

test95_GUMBEL <- VaRTest(alpha = 0.05, actual = V_test, VaR = VaR95_MC_GUMBEL)  
test99_GUMBEL <- VaRTest(alpha = 0.01, actual = V_test, VaR =VaR99_MC_GUMBEL)

test95_CLAYTON <- VaRTest(alpha = 0.05, actual = V_test, VaR = VaR95_MC_CLAYTON)  
test99_CLAYTON <- VaRTest(alpha = 0.01, actual = V_test, VaR =VaR99_MC_CLAYTON)

test95_FRANK <- VaRTest(alpha = 0.05, actual = V_test, VaR = VaR95_MC_FRANK)  
test99_FRANK <- VaRTest(alpha = 0.01, actual = V_test, VaR =VaR99_MC_FRANK)

test95_JOE <- VaRTest(alpha = 0.05, actual = V_test, VaR = VaR95_MC_JOE)  
test99_JOE <- VaRTest(alpha = 0.01, actual = V_test, VaR =VaR99_MC_JOE)


# ------------ Tables -----------------------------------------------------


# Loss Functions table

# t student
lopez95_T <- lopez_loss(V_test, VaR95_MC_T)
lopez99_T <- lopez_loss(V_test, VaR99_MC_T)

# Gauss
lopez99_GAUSS<- lopez_loss(V_test, VaR99_MC_GAUSS)
lopez95_GAUSS<- lopez_loss(V_test, VaR95_MC_GAUSS)

# Clayton
lopez99_CLAYTON<- lopez_loss(V_test, VaR99_MC_CLAYTON)
lopez95_CLAYTON <- lopez_loss(V_test, VaR95_MC_CLAYTON)

# FRANK
lopez99_FRANK<- lopez_loss(V_test, VaR99_MC_FRANK)
lopez95_FRANK <- lopez_loss(V_test, VaR95_MC_FRANK)

# GUMBEL
lopez99_GUMBEL<- lopez_loss(V_test, VaR99_MC_GUMBEL)
lopez95_GUMBEL <- lopez_loss(V_test, VaR95_MC_GUMBEL)

# JOE
lopez99_JOE<- lopez_loss(V_test, VaR99_MC_JOE)
lopez95_JOE <- lopez_loss(V_test, VaR95_MC_JOE)





# Copula names 
copula_names_final <- c("Gaussian", "t", "Clayton", "Frank", "Gumbel", "Joe")

#  Function to combine stat and p-value 
combine_stat_pval <- function(stat, pval) {
  ifelse(is.na(stat), "--", paste0(round(stat, 2), " (", round(pval, 4), ")"))
}

# === Confidence level 95% ===
ExpVio95 <- c(
  test95_GAUSS$expected.exceed,
  test95_T$expected.exceed,
  test95_CLAYTON$expected.exceed,
  test95_FRANK$expected.exceed,
  test95_GUMBEL$expected.exceed,
  test95_JOE$expected.exceed
)

ActVio95 <- c(
  test95_GAUSS$actual.exceed,
  test95_T$actual.exceed,
  test95_CLAYTON$actual.exceed,
  test95_FRANK$actual.exceed,
  test95_GUMBEL$actual.exceed,
  test95_JOE$actual.exceed
)

ZT95 <- ActVio95 / length(V_test)

LRuc_95_stat <- c(
  test95_GAUSS$uc.LRstat,
  test95_T$uc.LRstat,
  test95_CLAYTON$uc.LRstat,
  test95_FRANK$uc.LRstat,
  test95_GUMBEL$uc.LRstat,
  test95_JOE$uc.LRstat
)

LRuc_95_pval <- c(
  test95_GAUSS$uc.LRp,
  test95_T$uc.LRp,
  test95_CLAYTON$uc.LRp,
  test95_FRANK$uc.LRp,
  test95_GUMBEL$uc.LRp,
  test95_JOE$uc.LRp
)

LRcc_95_stat <- c(
  test95_GAUSS$cc.LRstat,
  test95_T$cc.LRstat,
  test95_CLAYTON$cc.LRstat,
  test95_FRANK$cc.LRstat,
  test95_GUMBEL$cc.LRstat,
  test95_JOE$cc.LRstat
)

LRcc_95_pval <- c(
  test95_GAUSS$cc.LRp,
  test95_T$cc.LRp,
  test95_CLAYTON$cc.LRp,
  test95_FRANK$cc.LRp,
  test95_GUMBEL$cc.LRp,
  test95_JOE$cc.LRp
)

Lopez95 <- c(
  lopez95_GAUSS,
  lopez95_T,
  lopez95_CLAYTON,
  lopez95_FRANK,
  lopez95_GUMBEL,
  lopez95_JOE
)

# Tabella finale 95%
tableVaR95 <- data.frame(
  Model = copula_names_final,
  ExpVio95 = ExpVio95,
  ActVio95 = ActVio95,
  MT95 = round(ZT95, 4),
  LRuc = combine_stat_pval(LRuc_95_stat, LRuc_95_pval),
  LRcc = combine_stat_pval(LRcc_95_stat, LRcc_95_pval),
  Lopez95 = round(Lopez95, 5)*(100/ActVio95)
)

# === Livello 99% ===
ExpVio99 <- c(
  test99_GAUSS$expected.exceed,
  test99_T$expected.exceed,
  test99_CLAYTON$expected.exceed,
  test99_FRANK$expected.exceed,
  test99_GUMBEL$expected.exceed,
  test99_JOE$expected.exceed
)

ActVio99 <- c(
  test99_GAUSS$actual.exceed,
  test99_T$actual.exceed,
  test99_CLAYTON$actual.exceed,
  test99_FRANK$actual.exceed,
  test99_GUMBEL$actual.exceed,
  test99_JOE$actual.exceed
)

ZT99 <- ActVio99 / length(V_test)

LRuc_99_stat <- c(
  test99_GAUSS$uc.LRstat,
  test99_T$uc.LRstat,
  test99_CLAYTON$uc.LRstat,
  test99_FRANK$uc.LRstat,
  test99_GUMBEL$uc.LRstat,
  test99_JOE$uc.LRstat
)

LRuc_99_pval <- c(
  test99_GAUSS$uc.LRp,
  test99_T$uc.LRp,
  test99_CLAYTON$uc.LRp,
  test99_FRANK$uc.LRp,
  test99_GUMBEL$uc.LRp,
  test99_JOE$uc.LRp
)

LRcc_99_stat <- c(
  test99_GAUSS$cc.LRstat,
  test99_T$cc.LRstat,
  test99_CLAYTON$cc.LRstat,
  test99_FRANK$cc.LRstat,
  test99_GUMBEL$cc.LRstat,
  test99_JOE$cc.LRstat
)

LRcc_99_pval <- c(
  test99_GAUSS$cc.LRp,
  test99_T$cc.LRp,
  test99_CLAYTON$cc.LRp,
  test99_FRANK$cc.LRp,
  test99_GUMBEL$cc.LRp,
  test99_JOE$cc.LRp
)

Lopez99 <- c(
  lopez99_GAUSS,
  lopez99_T,
  lopez99_CLAYTON,
  lopez99_FRANK,
  lopez99_GUMBEL,
  lopez99_JOE
)

# Tabella finale 99%
tableVaR99 <- data.frame(
  Model = copula_names_final,
  ExpVio99 = ExpVio99,
  ActVio99 = ActVio99,
  MT99 = round(ZT99, 4),
  LRuc = combine_stat_pval(LRuc_99_stat, LRuc_99_pval),
  LRcc = combine_stat_pval(LRcc_99_stat, LRcc_99_pval),
  Lopez99 = round(Lopez99, 5)*(100/ActVio99)
)















#----------------------------- Best copula Parameters stability ---------------------------------------

# Temporal variable
copula_param_df_t$time <- X_test$date[1:nrow(copula_param_df_t)]

series_cols <- names(copula_param_df_t)[names(copula_param_df_t) != "time"]

# Data for the correlation parameters
first_6_series <- copula_param_df_t %>%
  dplyr::select(time, all_of(series_cols[1:6])) %>%
  pivot_longer(cols = -time, names_to = "series", values_to = "value") %>%
  mutate(series = recode(series,
                         "rho.1" = "FTSEMIB-DAX",
                         "rho.2" = "FTSEMIB-CAC",
                         "rho.3" = "FTSEMIB-IBEX",
                         "rho.4" = "DAX-CAC",
                         "rho.5" = "DAX-IBEX",
                         "rho.6" = "CAC-IBEX"))

# Data for "degree of freedom"
last_series <- copula_param_df_t %>%
  dplyr::select(time, all_of(series_cols[length(series_cols)])) %>%
  pivot_longer(cols = -time, names_to = "series", values_to = "value")

# Chart Correlation coefficients
plot6_t <- ggplot(first_6_series, aes(x = time, y = value, color = series)) +
  geom_line(size = 0.8) +
  labs(title = "Correlation Coefficients Over Time (t-Copula)",
       x = "Time",
       y = "Value",
       color = "Pair") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text  = element_text(size = 20),             
    legend.title = element_text(size = 13, face = "bold"),
    plot.title   = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text    = element_text(size = 10),
    axis.title   = element_text(size = 12)
  ) +
  scale_color_brewer(type = "qual", palette = "Set1")

# Chart degree of freedom
plot_df_t <- ggplot(last_series, aes(x = time, y = value)) +
  geom_line(color = "darkblue", size = 0.8) +
  labs(title = "Degrees of Freedom (df) Over Time",
       x = "Time",
       y = "df") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text  = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

pdf("Data/Chart Post preliminar analysis/CopulaT_parameters.pdf", width = 10, height = 8)
grid.arrange(plot6_t, plot_df_t, ncol = 1, heights = c(2, 1))
dev.off()





#------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------- Parameters stability for marginals---------------------------------------------------------------------------------
extract_param <- function(param_name) {
  sapply(garch_param_df, function(x) {
    sapply(x, function(y) if (!is.null(y) && !is.na(y[param_name])) y[param_name] else NA)
  })
}

# Extract and transpose the matrices
df_matrix <- t(extract_param("shape"))
skew_matrix <- t(extract_param("skew"))

# Converst inot dataframe
df_shape <- as.data.frame(df_matrix)
df_skew <- as.data.frame(skew_matrix)
colnames(df_shape) <- c("ITA", "GER", "FRA", "SPA")
colnames(df_skew) <- c("ITA", "GER", "FRA", "SPA")
df_shape$time <- X_test$date[1:nrow(df_shape)]
df_skew$time <- X_test$date[1:nrow(df_skew)]

# Transform in long format and assing the proper labels
df_shape_long <- pivot_longer(df_shape, cols = -time, names_to = "Series", values_to = "Shape") %>%
  mutate(Series = recode(Series,
                         "ITA" = "FTSEMIB",
                         "GER" = "DAX40",
                         "FRA" = "CAC40",
                         "SPA" = "IBEX35"))

df_skew_long <- pivot_longer(df_skew, cols = -time, names_to = "Series", values_to = "Skewness") %>%
  mutate(Series = recode(Series,
                         "ITA" = "FTSEMIB",
                         "GER" = "DAX40",
                         "FRA" = "CAC40",
                         "SPA" = "IBEX35"))

# Chart Shape
p_shape <- ggplot(df_shape_long, aes(x = time, y = Shape, color = Series)) +
  geom_line(size = 0.8) +
  labs(title = "Time-varying 'Shape' (df) Parameter", x = "Date", y = "Shape") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

# Chart skewness
p_skew <- ggplot(df_skew_long, aes(x = time, y = Skewness, color = Series)) +
  geom_line(size = 0.8) +
  labs(title = "Time-varying 'Skewness' Parameter", x = "Date", y = "Skewness") +
  theme_minimal() +
  theme(
    legend.position = "none",  # ⬅️ rimuove la legenda
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

grid.arrange(p_shape, p_skew, ncol = 1)

pdf("Data/Chart Post preliminar analysis/Parameters_Marginals.pdf", width = 10, height = 8)
grid.arrange(p_shape, p_skew, ncol = 1)
dev.off()



coeff_variance <- data.frame(
  FTSEMIB_DAX  = sd(copula_param_df_t$rho.1, na.rm = TRUE) / mean(copula_param_df_t$rho.1, na.rm = TRUE),
  FTSEMIB_CAC  = sd(copula_param_df_t$rho.2, na.rm = TRUE) / mean(copula_param_df_t$rho.2, na.rm = TRUE),
  FTSEMIB_IBEX = sd(copula_param_df_t$rho.3, na.rm = TRUE) / mean(copula_param_df_t$rho.3, na.rm = TRUE),
  DAX_CAC      = sd(copula_param_df_t$rho.4, na.rm  = TRUE) / mean(copula_param_df_t$rho.4, na.rm = TRUE),
  DAX_IBEX     = sd(copula_param_df_t$rho.5, na.rm = TRUE) / mean(copula_param_df_t$rho.5, na.rm = TRUE),
  CAC_IBEX     = sd(copula_param_df_t$rho.6, na.rm = TRUE) / mean(copula_param_df_t$rho.6, na.rm = TRUE),
  
  shape_FTSEMIB = sd(df_shape$ITA, na.rm = TRUE) / mean(df_shape$ITA, na.rm = TRUE),
  shape_DAX     = sd(df_shape$GER, na.rm = TRUE) / mean(df_shape$GER, na.rm = TRUE),
  shape_CAC     = sd(df_shape$FRA, na.rm = TRUE) / mean(df_shape$FRA, na.rm = TRUE),
  shape_IBEX    = sd(df_shape$SPA, na.rm = TRUE) / mean(df_shape$SPA, na.rm = TRUE),
  
  skew_FTSEMIB  = sd(df_skew$ITA, na.rm = TRUE) / mean(df_skew$ITA, na.rm = TRUE),
  skew_DAX      = sd(df_skew$GER, na.rm = TRUE) / mean(df_skew$GER, na.rm = TRUE),
  skew_CAC      = sd(df_skew$FRA, na.rm = TRUE) / mean(df_skew$FRA, na.rm = TRUE),
  skew_IBEX     = sd(df_skew$SPA, na.rm = TRUE) / mean(df_skew$SPA, na.rm = TRUE)
)
coeff_variance <- t(coeff_variance*100)
colnames(coeff_variance) <- "CV"



