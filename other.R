# =============================================================================
#                           Other
# =============================================================================
# In this script we present:
# 1. News Impact Curve
# 2. Exceedence correlation (log returns and pseudo-observation)
# 3. Final tables (tables GARCH-Copula with Static Fit, GARCH-Copula with Dynamic Fit)
# 4. Chart ofo copula approach in parallel
# =============================================================================


# ========================== News Impact Curve ================================

series_names <- c("ITA", "GER", "FRA", "SPA")
series_labels <- c(ITA = "FTSEMIB", GER = "DAX", FRA = "CAC", SPA = "IBEX")

model_list <- list(
  sGARCH    = list(model = "sGARCH", submodel = NULL),
  AVGARCH   = list(model = "fGARCH", submodel = "AVGARCH"),
  gjrGARCH  = list(model = "gjrGARCH", submodel = NULL),
  TGARCH    = list(model = "fGARCH", submodel = "TGARCH")
)

model_colors <- c(sGARCH = "black", AVGARCH = "blue", gjrGARCH = "red", TGARCH = "darkgreen")
fit_results <- list()

# Models estimation
for (series in series_names) {
  fit_results[[series]] <- list()
  for (model_name in names(model_list)) {
    model_info <- model_list[[model_name]]
    
    if (model_name == "AVGARCH") {
      spec <- ugarchspec(
        variance.model = list(model = model_info$model,
                              garchOrder = c(1,1),
                              submodel = model_info$submodel),
        mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
        distribution.model = "sstd",
        fixed.pars = list(eta1 = 0, eta2 = 0)
      )
    } else {
      spec <- ugarchspec(
        variance.model = list(model = model_info$model,
                              garchOrder = c(1,1),
                              submodel = model_info$submodel),
        mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
        distribution.model = "sstd"
      )
    }
    
    fit <- ugarchfit(spec, data = X[[series]], solver = "hybrid")
    fit_results[[series]][[model_name]] <- fit
  }
}

# Chart
compute_nic <- function(model_fit, model_type) {
  p <- coef(model_fit)
  eps <- seq(-4, 4, length.out = 500)
  
  omega <- p["omega"]
  alpha <- p["alpha1"]
  beta  <- p["beta1"]
  
  if (model_type == "sGARCH") {
    nic <- sqrt(omega + alpha * eps^2 + beta)
  } else if (model_type == "AVGARCH") {
    nic <- omega + alpha * abs(eps) + beta
  } else if (model_type == "gjrGARCH") {
    gamma <- p["gamma1"]
    indicator <- ifelse(eps < 0, 1, 0)
    nic <- sqrt(omega + alpha * eps^2 + gamma * eps^2 * indicator + beta)
  } else if (model_type == "TGARCH") {
    eta1 <- p["eta11"]
    nic <- omega + alpha * (abs(eps) - eta1 * eps) + beta
  } else {
    stop("Unknown model type")
  }
  
  return(data.frame(eps = eps, nic = nic, model = model_type))
}

# Plot of a single serie
plot_nics_grid_ggplot <- function(series_name, fit_results, label_map) {
  plots <- list()
  
  for (model_name in names(fit_results[[series_name]])) {
    df <- compute_nic(fit_results[[series_name]][[model_name]], model_name)
    p <- ggplot(df, aes(x = eps, y = nic)) +
      geom_line(color = model_colors[[model_name]], linewidth = 1) +
      labs(title = paste(label_map[[series_name]], "-", model_name),
           x = expression(epsilon[t-1]), y = expression(sigma[t])) +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    plots[[model_name]] <- p
  }
  
  return(wrap_plots(plots, ncol = 2))
}

# All plots
for (series in series_names) {
  g <- plot_nics_grid_ggplot(series, fit_results, series_labels)
  ggsave(file.path("Data/Chart Post preliminar analysis", paste0("NIC_", series, ".pdf")),
         g, width = 8, height = 6)
  
}




# ====================== Excedeence correlation ===============================

U_matrix <- as.data.frame(U_matrix)
colnames(U_matrix) <- c("ITA", "GER", "FRA", "SPA")

# Function for bootstrap
plot_rho_matrix_gg_bootstrap <- function(data, title_prefix = "", B = 500) {
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  
  rho_extreme_boot <- function(X, Y, q, B = 500) {
    n <- length(X)
    qx <- quantile(X, q)
    qy <- quantile(Y, q)
    
    if (q <= 0.5) {
      idx <- which(X <= qx & Y <= qy)
    } else {
      idx <- which(X > qx & Y > qy)
    }
    
    if (length(idx) < 5) return(tibble(rho = NA, lower = NA, upper = NA, q = q))
    
    boots <- replicate(B, {
      samp_idx <- sample(idx, replace = TRUE)
      cor(X[samp_idx], Y[samp_idx], use = "complete.obs")
    })
    
    tibble(
      q = q,
      rho = mean(boots, na.rm = TRUE),
      lower = quantile(boots, 0.025, na.rm = TRUE),
      upper = quantile(boots, 0.975, na.rm = TRUE)
    )
  }
  
  index_labels <- c("ITA" = "FTSEMIB", "GER" = "DAX", "FRA" = "CAC", "SPA" = "IBEX")
  series_names <- c("ITA", "GER", "FRA", "SPA")
  pairs <- combn(series_names, 2, simplify = FALSE)
  q_vals <- seq(0.05, 0.95, by = 0.05)
  
  all_results <- map_dfr(pairs, function(pair) {
    X <- data[[pair[1]]]
    Y <- data[[pair[2]]]
    
    map_dfr(q_vals, function(q) {
      rho_extreme_boot(X, Y, q, B) %>%
        mutate(pair = paste(index_labels[pair[1]], "-", index_labels[pair[2]]))
    })
  })
  
  # Plot
  ggplot(all_results, aes(x = q, y = rho)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.4) +
    geom_line(color = "steelblue", linewidth = 0.8) +
    geom_point(color = "steelblue", size = 1.4) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    facet_wrap(~pair, ncol = 2) +
    labs(
      x = "Quantile (q)",
      y = expression(rho^e)
      # title = title_prefix  # rimosso titolo principale
    ) +
    theme_minimal(base_size = 13) +
    theme(
      # plot.title = element_text(size = 15, face = "bold", hjust = 0.5),  # disattivato
      strip.text = element_text(size = 30, face = "bold")
    )
}

p_log <- plot_rho_matrix_gg_bootstrap(X, "Log Returns")
ggsave("Data/Chart preliminar analysis/rho_extreme_logreturns_bootstrap.pdf", p_log, width = 12, height = 8)
p_pit <- plot_rho_matrix_gg_bootstrap(U_matrix, "PIT Residuals")
ggsave("Data/Chart Post preliminar analysis/rho_extreme_pit_bootstrap.pdf", p_pit, width = 12, height = 8)





# ========================= Final Tables ========================================


names(tableVaR95) <- c("Model", "ExpVio", "ActVio", "MT", "LRuc", "LRcc", "Lopez")
names(tableVaR99) <- c("Model", "ExpVio", "ActVio", "MT", "LRuc", "LRcc", "Lopez")

tableVaR_dynamic <- rbind(tableVaR95, tableVaR99)

names(tableVaR95_static) <- c("Model", "ExpVio", "ActVio", "MT", "LRuc", "LRcc", "Lopez")
names(tableVaR99_static) <- c("Model", "ExpVio", "ActVio", "MT", "LRuc", "LRcc", "Lopez")

tableVaR_static <- rbind(tableVaR95_static, tableVaR99_static)

names(tableVaR95_OT) <- c("Model", "ExpVio", "ActVio", "MT", "LRuc", "LRcc", "Lopez")
names(tableVaR99_OT) <- c("Model", "ExpVio", "ActVio", "MT", "LRuc", "LRcc", "Lopez")

tableVaR_OT <- rbind(tableVaR95_OT, tableVaR99_OT)


tableVaR95_all <- rbind(tableVaR95, tableVaR95_static, tableVaR95_OT)
tableVaR99_all <- rbind(tableVaR99, tableVaR99_static, tableVaR99_OT)



# ====================== PLOT Garch-Copula with static and dynamic fit ========


copula_names <- c("T", "GAUSS", "FRANK", "GUMBEL", "CLAYTON", "JOE")

# Separate lists 
plot_list_dynamic <- vector("list", length(copula_names))
plot_list_static  <- vector("list", length(copula_names))

# Create plot for dynamic fit
for (i in seq_along(copula_names)) {
  copula <- copula_names[i]
  VaR95 <- get(paste0("VaR95_MC_", copula))
  VaR99 <- get(paste0("VaR99_MC_", copula))
  
  p <- (function() {
    plot_VaR_hits(copula, V_test, X_test, VaR95, VaR99,
                  col = c("black", "darkblue", "lightblue", "red", "orange"))
    ggplot2::last_plot()
  })()
  
  p <- p + labs(title = paste(copula, "- Dynamic Fit"), x = NULL, y = NULL) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
          axis.text = element_text(size = 6))
  
  plot_list_dynamic[[i]] <- p
}

# Create plot for static fit
for (i in seq_along(copula_names)) {
  copula <- copula_names[i]
  VaR95 <- get(paste0("VaR95_MC_", copula, "_static"))
  VaR99 <- get(paste0("VaR99_MC_", copula, "_static"))
  
  p <- (function() {
    plot_VaR_hits(copula, V_test, X_test, VaR95, VaR99,
                  col = c("black", "darkgreen", "green", "red", "orange"))
    ggplot2::last_plot()
  })()
  
  p <- p + labs(title = paste(copula, "- Static Fit"), x = NULL, y = NULL) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
          axis.text = element_text(size = 6))
  
  plot_list_static[[i]] <- p
}

#  T, GAUSS, CLAYTON
indici_parte1 <- match(c("T", "GAUSS", "CLAYTON"), copula_names)
paired_plots_1 <- vector("list", length(indici_parte1) * 2)
for (i in seq_along(indici_parte1)) {
  idx <- indici_parte1[i]
  paired_plots_1[[2*i - 1]] <- plot_list_dynamic[[idx]]
  paired_plots_1[[2*i]]     <- plot_list_static[[idx]]
}
combined_part1 <- wrap_plots(paired_plots_1, ncol = 2, nrow = 3)

#  FRANK, GUMBEL, JOE
indici_parte2 <- match(c("FRANK", "GUMBEL", "JOE"), copula_names)
paired_plots_2 <- vector("list", length(indici_parte2) * 2)
for (i in seq_along(indici_parte2)) {
  idx <- indici_parte2[i]
  paired_plots_2[[2*i - 1]] <- plot_list_dynamic[[idx]]
  paired_plots_2[[2*i]]     <- plot_list_static[[idx]]
}
combined_part2 <- wrap_plots(paired_plots_2, ncol = 2, nrow = 3)

# Save
ggsave("Data/Chart Post preliminar analysis/COMPARISON_STATIC_vs_DYNAMIC_part1.pdf",
       combined_part1, width = 12, height = 9)

ggsave("Data/Chart Post preliminar analysis/COMPARISON_STATIC_vs_DYNAMIC_part2.pdf",
       combined_part2, width = 12, height = 9)
