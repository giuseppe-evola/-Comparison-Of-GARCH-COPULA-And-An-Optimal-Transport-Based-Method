# ==============================================================================
# Selection of Optimal ARMA-GARCH Models and A Priori Analysis of Copula Families
# ==============================================================================
# This script is structured in two main parts:
#
# 1. Selection of Best-Fitting ARMA-GARCH Models:
#    - For each stock index, we identify the optimal ARMA-GARCH specification
#      based on the function selection_optimal_arma_garch.
#
# 2. A Priori Analysis of Copula Families:
#    - Preliminary assessment of various static copula models to determine
#      which families are more suitable for capturing the dependence structure
#      observed in the standardized residuals of the selected marginal models.
# ==============================================================================


# ===================================================================== BEST GARCH MODELS =====================================================================


# Testing the best ARMA-GARCH model WITH MU0 COMPONENT
risultato_ITA1 <- selection_optimal_arma_garch(X$ITA)
risultato_GER1 <- selection_optimal_arma_garch(X$GER)
risultato_FRA1 <- selection_optimal_arma_garch(X$FRA)
risultato_SPA1 <- selection_optimal_arma_garch(X$SPA)

# It is easy to see that mu0 in not significant for all the time series analyzed. So we use re specify the models without this component.


# ITA [BEST MODEL]
spec_ITA <- ugarchspec(variance.model = list(model = "fGARCH", garchOrder = c(1, 1), submodel = "TGARCH"), # Specification
                       mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                       distribution.model = "sstd")
fit_ITA <- ugarchfit(spec_ITA, data = X$ITA, solver = "hybrid")                     # fitting

resid_standardized_ITA <- residuals(fit_ITA, standardize = TRUE)                  # standardized reisudals extraction
resid_standardized_ITA <- coredata(resid_standardized_ITA)

FinTS::ArchTest(resid_standardized_ITA, lags = 10)                                # LM-ARCH

Box.test(resid_standardized_ITA, lag = 20, type = c("Ljung-Box"))                # Ljung-box
acf(resid_standardized_ITA, main = "ACF of standardized residuals FTSEMIB")

Box.test(resid_standardized_ITA^2, lag = 20, type = c("Ljung-Box"))              # Ljung-Box ^2
acf(resid_standardized_ITA^2, main = "ACF of squared standardized residuals FTSEMIB")

pits_ITA <- rugarch::pdist("sstd", q = resid_standardized_ITA, shape = coef(fit_ITA)["shape"], skew = coef(fit_ITA)["skew"])
ks.test(pits_ITA, "punif")                                                       # Kolmogoroc Smirnov





# GER [BEST MODEL]
spec_GER <- ugarchspec(variance.model = list(model = "fGARCH", garchOrder = c(1, 1), submodel = "TGARCH"),
                       mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                       distribution.model = "sstd")
fit_GER <- ugarchfit(spec_GER, data = X$GER, solver = "hybrid")

resid_standardized_GER <- residuals(fit_GER, standardize = TRUE)
resid_standardized_GER <- coredata(resid_standardized_GER)

FinTS::ArchTest(resid_standardized_GER, lags = 10)

Box.test(resid_standardized_GER, lag = 20, type = c("Ljung-Box"))
acf(resid_standardized_GER, main = "ACF of standardized residuals DAX")

Box.test(resid_standardized_GER^2, lag = 20, type = c("Ljung-Box"))
acf(resid_standardized_GER^2, main = "ACF of squared standardized residuals DAX")

pits_GER <- rugarch::pdist("sstd", q = resid_standardized_GER, shape = coef(fit_GER)["shape"], skew = coef(fit_GER)["skew"])
ks.test(pits_GER, "punif")





# FRA [BEST MODEL]
spec_FRA <- ugarchspec(variance.model = list(model = "fGARCH", garchOrder = c(1, 1), submodel = "TGARCH"),
                       mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                       distribution.model = "sstd")
fit_FRA <- ugarchfit(spec_FRA, data = X$FRA, solver = "hybrid")

resid_standardized_FRA <- residuals(fit_FRA, standardize = TRUE)
FinTS::ArchTest(resid_standardized_FRA, lags = 10)

Box.test(resid_standardized_FRA, lag = 20, type = c("Ljung-Box"))
acf(resid_standardized_FRA, main = "ACF of standardized residuals CAC")

Box.test(resid_standardized_FRA^2, lag = 20, type = c("Ljung-Box"))
acf(resid_standardized_FRA^2, main = "ACF of squared standardized residuals CAC")

resid_standardized_FRA <- residuals(fit_FRA, standardize = TRUE)
resid_standardized_FRA <- coredata(resid_standardized_FRA)

pits_FRA <- rugarch::pdist("sstd", q = resid_standardized_FRA, shape = coef(fit_FRA)["shape"], skew = coef(fit_FRA)["skew"])
ks.test(pits_FRA, "punif")





# SPA [BEST MODEL]
spec_SPA <- ugarchspec(variance.model = list(model = "fGARCH", garchOrder = c(1, 1), submodel = "TGARCH"),
                      mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                      distribution.model = "sstd")
fit_SPA <- ugarchfit(spec_SPA, data = X$SPA, solver = "hybrid")

resid_standardized_SPA <- residuals(fit_SPA, standardize = TRUE)
resid_standardized_SPA <- coredata(resid_standardized_SPA)

FinTS::ArchTest(resid_standardized_SPA, lags = 10)

Box.test(resid_standardized_SPA, lag = 20, type = c("Ljung-Box"))
acf(resid_standardized_SPA, main = "ACF of standardized residuals CAC")

Box.test(resid_standardized_SPA^2, lag = 20, type = c("Ljung-Box"))
acf(resid_standardized_SPA^2, main = "ACF of squared standardized residuals CAC")

resid_standardized_SPA <- residuals(fit_SPA, standardize = TRUE)
resid_standardized_SPA <- coredata(resid_standardized_SPA)
pits_SPA <- rugarch::pdist("sstd", q = resid_standardized_SPA, shape = coef(fit_SPA)["shape"], skew = coef(fit_SPA)["skew"])
ks.test(pits_SPA, "punif")



#------------------------------------------------------------------------------------------------------
# Taylor Effect of log-returns (non standardized residuals)
acf(abs(residuals(fit_ITA, standardize = FALSE)), main = "ACF of ABS. non-standardized residuals FTSEMIB")
acf(abs(residuals(fit_GER, standardize = FALSE)), main = "ACF of ABS. non-standardized residuals DAX")
acf(abs(residuals(fit_FRA, standardize = FALSE)), main = "ACF of ABS. non-standardized residuals CAC")
acf(abs(residuals(fit_SPA, standardize = FALSE)), main = "ACF of ABS. non-standardized residuals IBEX")

acf(residuals(fit_ITA, standardize = FALSE)**2, main = "ACF of squared non-standardized residuals FTSEMIB")
acf(residuals(fit_GER, standardize = FALSE)**2, main = "ACF of squared non-standardized residuals DAX")
acf(residuals(fit_FRA, standardize = FALSE)**2, main = "ACF of squared non-standardized residuals CAC")
acf(residuals(fit_SPA, standardize = FALSE)**2, main = "ACF of squared non-standardized residuals IBEX")
#--------------------------------------------------------------------------------------------------

# QQ-Plots of standardized residuals

# FTSEMIB
qqplot_ITA <- ggplot() +
  stat_qq(aes(sample = pits_ITA), distribution = qunif, color = "blue") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1.2) +
  labs(title = "FTSEMIB", x = "Theoretical Quantiles", y = "Empirical Quantiles") +
  theme_minimal() +
  theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5))

# DAX
qqplot_GER <- ggplot() +
  stat_qq(aes(sample = pits_GER), distribution = qunif, color = "blue") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1.2) +
  labs(title = "DAX", x = "Theoretical Quantiles", y = "Empirical Quantiles") +
  theme_minimal() +
  theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5))

# CAC
qqplot_FRA <- ggplot() +
  stat_qq(aes(sample = pits_FRA), distribution = qunif, color = "blue") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1.2) +
  labs(title = "CAC", x = "Theoretical Quantiles", y = "Empirical Quantiles") +
  theme_minimal() +
  theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5))

# IBEX
qqplot_SPA <- ggplot() +
  stat_qq(aes(sample = pits_SPA), distribution = qunif, color = "blue") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1.2) +
  labs(title = "IBEX", x = "Theoretical Quantiles", y = "Empirical Quantiles") +
  theme_minimal() +
  theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5))

# Grid
grid.arrange(qqplot_ITA, qqplot_GER, qqplot_FRA, qqplot_SPA, ncol = 2)

ggsave(
  filename = "Data/Chart Post preliminar analysis/qqplot_pits_vs_uniform.pdf",
  plot = arrangeGrob(qqplot_ITA, qqplot_GER, qqplot_FRA, qqplot_SPA, ncol = 2),
  width = 10, height = 8
)




############# EXPORT STANDARDIZED RESIDUALS for Optima Transport phase #########
resids_df <- data.frame(
  FTSEMIB_resids = resid_standardized_ITA,
  DAX_resids     = resid_standardized_GER,
  CAC_resids     = resid_standardized_FRA,
  IBEX_resids    = resid_standardized_SPA
)
write.csv(resids_df, file = "standardized_residuals.csv", row.names = FALSE)


# Comparing ARMA-GARCH with mu0 vs ARMA-GARCH with no mu0

risultato_ITA1$table
risultato_GER1$table
risultato_FRA1$table
risultato_SPA1$table


risultato_ITA1$best_model
fit_ITA

risultato_GER1$best_model
fit_GER

risultato_FRA1$best_model
fit_FRA

risultato_SPA1$best_model
fit_SPA


##################################################################################################################################
############################################ table best specification for each time series #######################################
##################################################################################################################################

extract_parameters <- function(fit) {
  coefs <- round(coef(fit), 4)
  pvals <- round(fit@fit$matcoef[, 4], 4)
  parameters <- names(coefs)
  val <- sprintf("%.4f\n(%.4f)", coefs, pvals)
  names(val) <- parameters
  return(val)
}

# Extract parameter info for each index
params_ITA <- extract_parameters(fit_ITA)
params_GER <- extract_parameters(fit_GER)
params_FRA <- extract_parameters(fit_FRA)
params_SPA <- extract_parameters(fit_SPA)

# Get the union of all parameter names (some models may have different components)
all_parameters <- union(union(names(params_ITA), names(params_GER)),
                        union(names(params_FRA), names(params_SPA)))

# Helper function to ensure all parameters are present for each model
align_parameters <- function(params, all_names) {
  sapply(all_names, function(p) ifelse(p %in% names(params), params[p], ""))
}

# Build the summary table
param_table <- tibble(
  Parameter = all_parameters,
  FTSEMIB   = align_parameters(params_ITA, all_parameters),
  DAX       = align_parameters(params_GER, all_parameters),
  CAC40     = align_parameters(params_FRA, all_parameters),
  IBEX35    = align_parameters(params_SPA, all_parameters)
)

# Retrieve the name of the fitted model for each index
model_name <- function(fit) {
  spec <- fit@model$modeldesc
  arma <- paste0("ARMA(", paste(fit@model$arma[1:2], collapse = ","), ")")
  paste(arma, "-", spec$vmodel, "(1,1)", "-", spec$distribution)
}

# Custom column headers including the model specification
colnames(param_table)[2:5] <- c(
  paste0("FTSEMIB\n", model_name(fit_ITA)),
  paste0("DAX\n", model_name(fit_GER)),
  paste0("CAC40\n", model_name(fit_FRA)),
  paste0("IBEX35\n", model_name(fit_SPA))
)

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################



# ===================================================================== A-PRIORI ANALYSIS - COPULA MODELS =====================================================================

# Df of Pseudo observation
U_matrix <- cbind(ITA = pits_ITA,
                  GER = pits_GER,
                  FRA = pits_FRA,
                  SPA = pits_SPA)



# Copula definition [Estimation method: IMF]
fit.gaussCopula <- fitCopula(normalCopula(dim = 4, dispstr = "un"),  data = U_matrix, method = "ml")
fit.tCopula <- fitCopula(tCopula(dim = 4, dispstr = "un"), data = U_matrix, method = "ml") 
fit.claytonCopula <- fitCopula(claytonCopula(dim = 4), data = U_matrix, method = "ml")
fit.frankCopula <- fitCopula(frankCopula(dim = 4), data = U_matrix, method = "ml")
fit.gumbelCopula <- fitCopula(gumbelCopula(dim = 4), data = U_matrix, method = "ml")
fit.joeCopula <- fitCopula(joeCopula(dim=4), data = U_matrix, method = "ml")

lambda_gauss <- lambda(fit.gaussCopula@copula)
lambda_t     <- lambda(fit.tCopula@copula)
lambda_clayton <- lambda(fit.claytonCopula@copula)
lambda_frank     <- lambda(fit.frankCopula@copula)
lambda_gumbel <- lambda(fit.gumbelCopula@copula)
lambda_joe    <- lambda(fit.joeCopula@copula)





# ---- Table of Elliptical and Archimedean Copulas ----
elliptical_table <- bind_rows(
  get_elliptical_info(fit.gaussCopula, "gaussian"),
  get_elliptical_info(fit.tCopula, "t")
) %>% relocate(Copula)

lambda(fit.gaussCopula@copula)

lambda_table <- data.frame(
  Gauss = lambda(fit.gaussCopula@copula),
  t     = lambda(fit.tCopula@copula)
)

rownames(lambda_table) <- c("FTSEMIB-DAX - Low", "FTSEMIB-CAC - Low", "FTSEMIB-IBEX - Low",
                 "DAX-CAC - Low", "DAX-IBEX - Low", "CAC-IBEX - Low",
                 "FTSEMIB-DAX - Up", "FTSEMIB-CAC - Up", "FTSEMIB-IBEX - Up",
                 "DAX-CAC - Up", "DAX-IBEX - Up", "CAC-IBEX - Up")

archimedean_table <- bind_rows(
  get_archimedean_info(fit.claytonCopula, "clayton"),
  get_archimedean_info(fit.frankCopula, "frank"),
  get_archimedean_info(fit.gumbelCopula, "gumbel"),
  get_archimedean_info(fit.joeCopula, "joe")  # 
) %>% relocate(Copula)

rownames(elliptical_table) <- NULL  
rownames(archimedean_table) <- NULL 

lambda_table_arch <- data.frame(
  Clayton    = c(lambda(fit.claytonCopula@copula)[1], lambda(fit.claytonCopula@copula)[2]),
  Frank  = c(lambda(fit.frankCopula@copula)[1], lambda(fit.frankCopula@copula)[2]),
  Gumbel = c(lambda(fit.gumbelCopula@copula)[1], lambda(fit.gumbelCopula@copula)[2]),
  Joe =  c(lambda(fit.joeCopula@copula)[1], lambda(fit.joeCopula@copula)[2])
)




# ------------------------------- VISUAL INSPECTION for BEST COPULA ----------------------------------------


# Pseudo-observation
df_pit <- as.data.frame(U_matrix)
colnames(df_pit) <- c("FTSEMIB", "DAX", "CAC", "IBEX")

# Create the couples combination
pair_names <- combn(colnames(df_pit), 2, simplify = FALSE)

# Create a panel for each couple
make_plot <- function(var1, var2, data) {
  tau_val <- round(cor(data[[var1]], data[[var2]], method = "kendall"), 3)
  
  ggplot(data, aes_string(x = var1, y = var2)) +
    geom_point(alpha = 0.3, size = 0.7) +
    annotate("text", x = 0.05, y = 0.95, 
             label = paste0("τ = ", tau_val), 
             hjust = 0, vjust = 1, color = "firebrick", size = 10, fontface = "bold") +
    theme_minimal(base_size = 11) +
    theme(axis.title = element_text(size = 20))+
    coord_fixed() +
    labs(x = var1, y = var2)
}

# Create 6 plots
plot_list <- lapply(pair_names, function(pair) make_plot(pair[1], pair[2], df_pit))

#Grid 2x3
final_plot <- gridExtra::grid.arrange(grobs = plot_list, ncol = 3)

ggsave(
  filename = "Data/Chart Post preliminar analysis/scatter_pseudo_observations_kendall_tau.pdf",
  plot = gridExtra::arrangeGrob(grobs = plot_list, ncol = 3),
  width = 12, height = 8
)






# SIMULATIONS FROM FITTED COPULAS vs REAL MARGINS (STANDARDIZED RESIDUALS PITs)
index_labels <- c("ITA" = "FTSEMIB", "GER" = "DAX", "FRA" = "CAC", "SPA" = "IBEX")

# Fitted copulas
copula_list <- list(
  Gaussian = fit.gaussCopula@copula,
  t        = fit.tCopula@copula,
  Clayton  = fit.claytonCopula@copula,
  Frank    = fit.frankCopula@copula,
  Gumbel   = fit.gumbelCopula@copula,
  Joe      = fit.joeCopula@copula
)

series_names <- colnames(U_matrix)
pair_names <- combn(series_names, 2, simplify = FALSE)
pair_labels <- sapply(pair_names, function(pair) paste(index_labels[pair[1]], "vs", index_labels[pair[2]]))

# matrix creation
all_plot_matrix <- matrix(list(), nrow = length(copula_list), ncol = length(pair_names))

# Fill the matrix
copula_names <- names(copula_list)
for (r in seq_along(copula_names)) {
  cop_name <- copula_names[r]
  cop <- copula_list[[cop_name]]
  sim_data <- rCopula(1000, cop)
  
  for (c in seq_along(pair_names)) {
    pair <- pair_names[[c]]
    var1 <- pair[1]
    var2 <- pair[2]
    
    idx1 <- which(series_names == var1)
    idx2 <- which(series_names == var2)
    
    df_real <- data.frame(U1 = U_matrix[, var1], U2 = U_matrix[, var2])
    df_sim  <- data.frame(U1 = sim_data[, idx1], U2 = sim_data[, idx2])
    
    p <- ggplot() +
      geom_point(data = df_real, aes(x = U1, y = U2),
                 color = "blue", alpha = 0.4, size = 1.2) +
      geom_point(data = df_sim, aes(x = U1, y = U2),
                 color = "#ff7f0e", alpha = 0.5, size = 0.9) +
      theme_minimal(base_size = 9) +
      theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        plot.margin = margin(2, 2, 2, 2)
      )
    
    all_plot_matrix[[r, c]] <- p
  }
}

# Converts in list for grid.arrange
plot_list <- as.vector(t(all_plot_matrix))

# Labels
col_title_grobs <- lapply(pair_labels, function(lab) {
  textGrob(lab, gp = gpar(fontsize = 12, fontface = "bold"))
})

row_title_grobs <- lapply(names(copula_list), function(lab) {
  textGrob(lab, rot = 90, gp = gpar(fontsize = 12, fontface = "bold"))
})

# Save
pdf("Data/Chart Post preliminar analysis/panels_copula_simulations_matrix.pdf", width = 18, height = 16)

# Create the wanted layout
pushViewport(viewport(layout = grid.layout(7, 7, 
                                           heights = unit(c(2, rep(1, 6)), c("cm", rep("null", 6))),
                                           widths = unit(c(2.5, rep(1, 6)), c("cm", rep("null", 6))))))

# Column Labels
for(i in 1:6) {
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = i+1))
  grid.text(pair_labels[i], gp = gpar(fontsize = 20, fontface = "bold"))
  popViewport()
}

# Rows Labels
for(i in 1:6) {
  pushViewport(viewport(layout.pos.row = i+1, layout.pos.col = 1))
  grid.text(names(copula_list)[i], rot = 90, gp = gpar(fontsize = 20, fontface = "bold"))
  popViewport()
}

# Add the charts
for(i in 1:6) {
  for(j in 1:6) {
    pushViewport(viewport(layout.pos.row = i+1, layout.pos.col = j+1))
    plot_idx <- (i-1)*6 + j
    print(plot_list[[plot_idx]], newpage = FALSE)
    popViewport()
  }
}

popViewport()
dev.off()
