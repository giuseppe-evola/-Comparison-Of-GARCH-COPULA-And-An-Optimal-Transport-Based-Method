# =====================================================================
# Functions - Master thesis
# =====================================================================
# In this script you an find some of the functions used in this work
# =====================================================================


# ---------Select best ARMA-GARCH model -----------
selection_optimal_arma_garch <- function(series, alpha = .05) {
  if(!is.ts(series) && !is.numeric(series)) stop("Input must be a time series or numeric vector")
  if(!is.ts(series)) series <- as.ts(series)
  
  models <- list(
    list(type = "sGARCH", submodel = NULL),
    list(type = "gjrGARCH", submodel = NULL),
    list(type = "fGARCH", submodel = "TGARCH"),
    list(type = "fGARCH", submodel = "AVGARCH")
  )
  
  distributions <- c("norm", "std", "sstd")
  p_orders <- 0:2
  q_orders <- 0:2
  
  all_models <- list()
  model_count <- 0
  
  for(model_info in models) {
    for(p in p_orders) {
      for(q in q_orders) {
        for(dist_type in distributions) {
          model_type <- model_info$type
          submodel <- model_info$submodel
          model_name <- ifelse(
            is.null(submodel),
            paste0(model_type, "-", dist_type, "-ARMA(", p, ",", q, ")"),
            paste0(model_type, "-", submodel, "-", dist_type, "-ARMA(", p, ",", q, ")")
          )
          
          cat("\nEstimating model:", model_name, "\n")
          
          tryCatch({
            if (!is.null(submodel) && model_type == "fGARCH" && submodel == "AVGARCH") {
              spec <- ugarchspec(
                variance.model = list(model = model_type, garchOrder = c(1, 1), submodel = submodel),
                mean.model = list(armaOrder = c(p, q), include.mean = TRUE),
                distribution.model = dist_type,
                fixed.pars = list(eta11 = 0, eta21 = 0)
              )
            } else {
              spec <- ugarchspec(
                variance.model = list(model = model_type, garchOrder = c(1, 1), submodel = submodel),
                mean.model = list(armaOrder = c(p, q), include.mean = TRUE),
                distribution.model = dist_type
              )
            }
            
            cat("Specification created. Starting estimation...\n")
            fit <- ugarchfit(spec = spec, data = series, solver = "hybrid")
            
            model_count <- model_count + 1
            bic_value <- infocriteria(fit)[2, 1]
            
            all_models[[model_count]] <- list(
              fit = fit,
              name = model_name,
              bic = bic_value,
              distribution = dist_type
            )
            
            cat("BIC:", bic_value, "\n")
          }, error = function(e) {
            cat("Estimation error:", e$message, "\n")
          })
        }
      }
    }
  }
  
  if(model_count == 0) {
    cat("No model was successfully estimated!\n")
    return(NULL)
  }
  
  cat("\nSuccessfully estimated models:", model_count, "out of", length(models)*length(p_orders)*length(q_orders)*length(distributions), "models tested.\n")
  
  # Ordina per BIC crescente
  bic_values <- sapply(all_models, function(x) x$bic)
  order_idx <- order(bic_values)
  sorted_models <- all_models[order_idx]
  top_n <- min(4, length(sorted_models))
  top_models <- sorted_models[1:top_n]
  
  summary_table <- data.frame(Model = character(),
                              AIC = numeric(),
                              BIC = numeric(),
                              ARCH_LM = character(),
                              Ljung_Box = character(),
                              KS_Test = character(),
                              stringsAsFactors = FALSE)
  
  best_model <- NULL
  valid_model_found <- FALSE
  
  for (model in top_models) {
    fit <- model$fit
    std_resid <- residuals(fit, standardize = TRUE)
    
    aic_val <- infocriteria(fit)[1, 1]
    bic_val <- round(model$bic, 4)
    
    arch_test <- FinTS::ArchTest(std_resid, lags = 10)
    arch_stat <- round(arch_test$statistic, 4)
    arch_pval <- round(arch_test$p.value, 4)
    
    lb_test <- Box.test(std_resid, lag = 20, type = "Ljung-Box")
    lb_stat <- round(lb_test$statistic, 4)
    lb_pval <- round(lb_test$p.value, 4)
    
    dist <- model$distribution
    if (dist == "norm") {
      pit_values <- pnorm(std_resid)
    } else if (dist == "std") {
      shape <- coef(fit)["shape"]
      pit_values <- pt(std_resid, df = shape)
    } else if (dist == "sstd") {
      shape <- coef(fit)["shape"]
      skew <- coef(fit)["skew"]
      pit_values <- rugarch::pdist("sstd", q = std_resid, shape = shape, skew = skew)
    }
    
    ks_test <- ks.test(pit_values, "punif")
    ks_stat <- round(ks_test$statistic, 4)
    ks_pval <- round(ks_test$p.value, 4)
    
    summary_table <- rbind(summary_table, data.frame(
      Model = model$name,
      AIC = round(aic_val, 4),
      BIC = bic_val,
      ARCH_LM = paste0(arch_stat, " (", arch_pval, ")"),
      Ljung_Box = paste0(lb_stat, " (", lb_pval, ")"),
      KS_Test = paste0(ks_stat, " (", ks_pval, ")"),
      stringsAsFactors = FALSE
    ))
    
    if (!valid_model_found && all(c(arch_pval, lb_pval, ks_pval) >= alpha)) {
      best_model <- fit
      valid_model_found <- TRUE
    }
  }
  
  # rank by bic
  summary_table <- summary_table[order(summary_table$BIC), ]
  
  return(list(
    best_model = best_model,
    table = summary_table
  ))
}

risultato_ITA1 <- selection_optimal_arma_garch(X$ITA)
risultato_GER1 <- selection_optimal_arma_garch(X$GER)
risultato_FRA1 <- selection_optimal_arma_garch(X$FRA)
risultato_SPA1 <- selection_optimal_arma_garch(X$SPA)

#---------------------------------------------------------------
# Function to transform standardized residuals into a dataframe with uniform distribution
uniform_trasformer <- function(fits) {
  
  pit_list <- list()
  
  for (i in seq_along(fits)) {
    fit <- fits[[i]]
    
    # Check the class (case insensitive)
    if (!any(tolower(class(fit)) == "ugarchfit")) {
      stop(paste("The object at position", i, "is not a valid uGARCHfit object."))
    }
    
    # Extract standardized residuals
    res <- residuals(fit, standardize = TRUE)
    
    # Basic check: residuals should not be constant or NA
    if (any(is.na(res)) || sd(res, na.rm = TRUE) == 0) {
      stop(paste("Invalid residuals at position", i))
    }
    
    # Get distribution information
    dist_name <- fit@model$modeldesc$distribution
    coef <- fit@fit$coef
    shape <- ifelse("shape" %in% names(coef), coef["shape"], NA)
    skew  <- ifelse("skew"  %in% names(coef),  coef["skew"],  NA)
    
    # Apply PIT transformation using the same distribution as in the model
    # For standardized residuals
    if (dist_name == "norm") {
      u <- pdist(q = res, dist = "norm")
    } else if (dist_name == "std") {
      u <- pdist(q = res, dist = "std", shape = shape)
    } else if (dist_name == "sstd") {
      u <- pdist(q = res, dist = "sstd", shape = shape, skew = skew)
    } else {
      stop(paste("Unsupported distribution:", dist_name))
    }
    
    # Clipping extreme values to avoid 0 or 1 (numerical issues in qdist later)
    u <- pmin(pmax(u, 1e-8), 1 - 1e-8)
    
    pit_list[[i]] <- u
  }
  
  # Create dataframe with dynamic column names
  pit_df <- as.data.frame(pit_list)
  colnames(pit_df) <- paste0("PIT", seq_along(fits))
  
  return(pit_df)
}



#----------------------------------------------------------------------------

# Function to compute the simulated log-returns after having the simulated PITs
simulate_logreturns <- function(fit_garch, zt_uniform_sim) {
  # Get the forecast for the next period (t+1)
  forecast <- ugarchforecast(fit_garch, n.ahead = 1)
  vol_forecast <- as.numeric(sigma(forecast))
  mean_forecast <- as.numeric(fitted(forecast))  # AR(p) and mu components are included here, if present
  
  zt_uniform_sim <- as.vector(zt_uniform_sim)
  
  # Ensure PITs are strictly within (0,1) to avoid Inf in qdist
  zt_uniform_sim <- pmin(pmax(zt_uniform_sim, 1e-8), 1 - 1e-8)
  
  # Get the distribution parameters
  if("shape" %in% names(coef(fit_garch))) {
    shape_param <- coef(fit_garch)["shape"]
  } else {
    shape_param <- NA
  }
  
  if("skew" %in% names(coef(fit_garch))) {
    skew_param <- coef(fit_garch)["skew"]
  } else {
    skew_param <- NA
  }
  
  # Get the distribution type
  dist_type <- fit_garch@model$modeldesc$distribution
  
  # Transform uniform PITs to distribution-specific innovations using qdist
  if (dist_type == "norm") {
    zt <- qdist("norm", p = zt_uniform_sim, mu = 0, sigma = 1)
  } else if (dist_type == "std") {
    zt <- qdist("std", p = zt_uniform_sim, mu = 0, sigma = 1, shape = shape_param)
  } else if (dist_type == "sstd") {
    zt <- qdist("sstd", p = zt_uniform_sim, mu = 0, sigma = 1, shape = shape_param, skew = skew_param)
  } else {
    stop(paste("Unsupported distribution:", dist_type))
  }
  
  # Compute the simulated return
  MC_logrets <- mean_forecast + (zt * vol_forecast)
  
  return(MC_logrets)
}


# ---------------------------------------------------------

# LOPEZ LOSS FUNCTIONS

lopez_loss <- function(PL, VaR) {
  res <- sum(ifelse(PL < VaR,  abs(VaR - PL), 0))
  return(res)
}


#-------------------------------------------------------------------
# CHARTS for each copula
plot_VaR_hits <- function(copula_name, V_test, X_test, VaR95, VaR99, col = c("black", "blue", "lightblue","red", "orange")) {
  df <- data.frame(
    Date = X_test[["date"]],
    Return = -V_test,
    VaR95 = -VaR95,
    VaR99 = -VaR99
  )
  
  # Violations
  df$Hit_99 <- df$Return > df$VaR99
  df$Hit_95 <- df$Return > df$VaR95
  
  # Points
  df$Col <- NA
  df$Col[df$Hit_99] <- col[4]               # Red per VaR99
  df$Col[!df$Hit_99 & df$Hit_95] <- col[5] # Yellow per solo VaR95
  
  p <- ggplot(df, aes(x = Date)) +
    geom_line(aes(y = Return), color = col[1], size = 0.7) +
    geom_line(aes(y = VaR95), color = col[3], size = 0.7) +
    geom_line(aes(y = VaR99), color = col[2], size = 0.7) +
    
    #VaR95
    geom_point(data = subset(df, Col == col[5]), aes(y = Return),
               color = col[5], size = 2.5) +
    
    # VaR99
    geom_point(data = subset(df, Col == col[4]), aes(y = Return),
               color = col[4], size = 3.5) +
    
    labs(title = paste("VaR 95% and 99% - Copula", copula_name),
         y = "Portfolio Return",
         x = "Date") +
    theme_bw(base_size = 10)
  
  print(p)
}

# VaR OT
plot_VaR_OT <- function(V_test, dates, VaR95, VaR99, 
                        col = c("black", "purple", "green", "red", "orange"), 
                        title_text = "VaR 95% and 99% - Optimal Transport") {
  
  df <- data.frame(
    Date = dates,
    Return = -V_test,     # Loss
    VaR95 = -VaR95,
    VaR99 = -VaR99
  )
  
  # Violazioni (loss > VaR)
  df$Hit_99 <- df$Return > df$VaR99
  df$Hit_95 <- df$Return > df$VaR95
  
  # Colorazione
  df$Col <- NA
  df$Col[df$Hit_99] <- col[4]              # rosso
  df$Col[!df$Hit_99 & df$Hit_95] <- col[5] # arancio
  
  # Plot aggiornato
  p <- ggplot(df, aes(x = Date)) +
    geom_line(aes(y = Return), color = col[1], size = 0.7) +
    geom_line(aes(y = VaR95),   color = col[3], size = 0.7) +
    geom_line(aes(y = VaR99),   color = col[2], size = 0.7) +
    geom_point(data = subset(df, Col == col[5]), aes(y = Return),
               color = col[5], size = 2.5) +
    geom_point(data = subset(df, Col == col[4]), aes(y = Return),
               color = col[4], size = 3.5) +
    labs(title = title_text,
         y = NULL,         # <-- rimuove etichetta asse y
         x = "Date") +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      axis.title.y = element_blank()
    )
  
  print(p)
}
#----------------------------------


# Tail dependence functions
#lower_tail_dep <- function(u, v, q = 0.05) {
#  mean((u < q) & (v < q), na.rm = TRUE) / q
#}

#upper_tail_dep <- function(u, v, q = 0.95) {
#  mean((u > q) & (v > q), na.rm = TRUE) / (1 - q)
#}

#-----------------------------------------------

# Copulas tables
get_elliptical_info <- function(fit, name) {
  tibble(
    Copula = name,
    LogLik = logLik(fit)[1],
    AIC    = AIC(fit),
    BIC    = BIC(fit),
    Params = paste(round(coef(fit), 4), collapse = ", ")
  )
}

# Copulas tables
get_archimedean_info <- function(fit, name) {
  tibble(
    Copula = name,
    LogLik = logLik(fit)[1],
    AIC    = AIC(fit),
    BIC    = BIC(fit),
    Param  = round(coef(fit), 4)
  )
}


