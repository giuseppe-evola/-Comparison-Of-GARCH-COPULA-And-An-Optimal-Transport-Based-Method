# =====================================================================
# Preliminary Analysis Prior to Data Modelling - Master thesis
# =====================================================================
# This script performs the following steps:
# 1. Split the dataset into training and test sets
# 2. Compute and tabulate descriptive statistics for each asset
# 3. Explore and visualize stylized facts observed in financial time series,
#    as discussed in the thesis (e.g., volatility clustering, heavy tails, etc.)
# =====================================================================


#Libraries (all imported in this script)
library(dplyr)
library(ggplot2)
library(xts)
library(tidyr)
library(moments)
library(tseries)
library(forecast)
library(MASS)
library(copula)
library(rugarch)
library(xts)
library(nortest)
library(fitdistrplus)
library(stats)
library(tseries)
library(FinTS)
library(GGally)
library(gridExtra)
library(profvis)
library(parallel)
library(patchwork)
library(zoo)
library(tidyverse)
library(patchwork)
library(ggcorrplot)
library(grid)
library(knitr)
library(kableExtra)
setwd("C:/Users/giuse/OneDrive/Desktop/Preparazione tesi/Thesis Code")


# Data manipulation

italy   <- read.csv("Data/italy.csv", stringsAsFactors = FALSE)
germany <- read.csv("Data/germany.csv", stringsAsFactors = FALSE)
france  <- read.csv("Data/france.csv", stringsAsFactors = FALSE)
spain   <- read.csv("Data/spain.csv", stringsAsFactors = FALSE)

italy <- italy[1:nrow(italy),1:2]
colnames(italy) <- c("date", "price_italy")
italy$date <- as.Date(italy$date, format = "%m/%d/%Y")
italy$price_italy <- as.numeric(gsub(",", "", italy$price_italy))

germany <- germany[1:nrow(germany),1:2]
colnames(germany) <- c("date", "price_germany")
germany$date <- as.Date(germany$date, format = "%m/%d/%Y")
germany$price_germany <- as.numeric(gsub(",", "", germany$price_germany))

france <- france[1:nrow(france),1:2]
colnames(france) <- c("date", "price_france")
france$date <- as.Date(france$date, format = "%m/%d/%Y")
france$price_france <- as.numeric(gsub(",", "", france$price_france))

spain <- spain[1:nrow(spain),1:2]
colnames(spain) <- c("date", "price_spain")
spain$date <- as.Date(spain$date, format = "%m/%d/%Y")
spain$price_spain <- as.numeric(gsub(",", "", spain$price_spain))

italy$date <- as.Date(italy$date, format = "%d-%m-%y")
germany$date   <- as.Date(germany$date, format = "%d-%m-%y")
france$date    <- as.Date(france$date, format = "%d-%m-%y")
spain$date      <- as.Date(spain$date, format = "%d-%m-%y")

merged <- merge(italy, germany, by = "date")
merged <- merge(merged, france, by = "date")
df <- merge(merged, spain, by = "date")


str(df)

# Using Log-rets
X_tot <- data.frame(date = df$date[-1],  # remoe first date
                         ITA = diff(log(df$price_italy)),
                         GER = diff(log(df$price_germany)),
                        FRA= diff(log(df$price_france)),
                         SPA = diff(log(df$price_spain)))
str(X_tot)
X_tot[nrow(X_tot),]


# We divide into train set and test set:
# Train test: first  approximaltely 6 years of data (1500 datapoints)
# Test set: the other 1110 datapoints

# Test Set
X_test <- X_tot[1501:nrow(X_tot),]
nrow(X_test)

# Train set
X <- X_tot[1:1500,]
nrow(X)



#--------------------------------------------------- Basic Statistics ------------------------------------------------------------------------


# Descriptive statistics
stats_table <- data.frame(
  Statistic = c("Min", "FirstQuart", "Median", "Mean", "ThirdQuart", "Max", "SD", "Skewness", "ExcKurtosis"),
  ITA = c(min(X$ITA),quantile(X$ITA, 0.25), median(X$ITA), mean(X$ITA), quantile(X$ITA, 0.75), max(X$ITA), sd(X$ITA), skewness(X$ITA), kurtosis(X$ITA)),
  GER = c(min(X$GER),quantile(X$GER, 0.25), median(X$GER), mean(X$GER), quantile(X$GER, 0.75), max(X$GER), sd(X$GER), skewness(X$GER), kurtosis(X$GER)),
  FRA = c(min(X$FRA),quantile(X$FRA, 0.25), median(X$FRA), mean(X$FRA), quantile(X$FRA, 0.75), max(X$FRA), sd(X$FRA), skewness(X$FRA), kurtosis(X$FRA)),
  SPA = c(min(X$SPA),quantile(X$SPA, 0.25), median(X$SPA), mean(X$SPA), quantile(X$SPA, 0.75), max(X$SPA), sd(X$SPA), skewness(X$SPA), kurtosis(X$SPA))
)

stats_table[, 2:5] <- round(stats_table[, 2:5], 2)

# Functions for the p values formatting
fmt_test <- function(test_result) {
  stat <- round(test_result$statistic, 2)
  pval <- format.pval(test_result$p.value, digits = 3, eps = .001)
  paste0(stat, " (", pval, ")")
}

# ADD tests
tests <- data.frame(
  Statistic = c("JB", "LB", "LB$^2$", "ADF"),
  ITA = c(fmt_test(jarque.bera.test(X$ITA)),
          fmt_test(Box.test(X$ITA, type = "Ljung-Box", lag = 20)),
          fmt_test(Box.test(X$ITA^2, type = "Ljung-Box", lag = 20)),
          fmt_test(adf.test(X$ITA))),
  GER = c(fmt_test(jarque.bera.test(X$GER)),
          fmt_test(Box.test(X$GER, type = "Ljung-Box", lag = 20)),
          fmt_test(Box.test(X$GER^2, type = "Ljung-Box", lag = 20)),
          fmt_test(adf.test(X$GER))),
  FRA = c(fmt_test(jarque.bera.test(X$FRA)),
          fmt_test(Box.test(X$FRA, type = "Ljung-Box", lag = 20)),
          fmt_test(Box.test(X$FRA^2, type = "Ljung-Box", lag = 20)),
          fmt_test(adf.test(X$FRA))),
  SPA = c(fmt_test(jarque.bera.test(X$SPA)),
          fmt_test(Box.test(X$SPA, type = "Ljung-Box", lag = 20)),
          fmt_test(Box.test(X$SPA^2, type = "Ljung-Box", lag = 20)),
          fmt_test(adf.test(X$SPA)))
)

# Combine to the final table
final_table <- rbind(stats_table, tests)




# --------------------------------------------- Normality Check ----------------------------------------------------------------

# Path
base_path <- "Data/Chart preliminar analysis"

index_colors <- c(
  "GER" = "black",
  "SPA" = "red",
  "FRA" = "blue",
  "ITA" = "green"
)

# Log-rets in "long format"
X_long <- X %>%
  pivot_longer(-date, names_to = "Index", values_to = "Log_Return")

# Price index
price_index <- X %>%
  mutate(across(-date, ~ cumprod(1 + .), .names = "{.col}")) %>%
  pivot_longer(-date, names_to = "Index", values_to = "Price_Index")

# Couple the charts
plot_pair <- function(index_name, index_label) {
  # Log-Returns
  p1 <- ggplot(filter(X_long, Index == index_name), aes(x = date, y = Log_Return)) +
    geom_line(color = index_colors[index_name], size = 0.8) +
    labs(title = paste("Log-Returns -", index_label),
         x = "Date", y = "Log-Returns") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(size = 30, face = "bold"))
  
  # Price Index
  p2 <- ggplot(filter(price_index, Index == index_name), aes(x = date, y = Price_Index)) +
    geom_line(color = index_colors[index_name], size = 0.8) +
    labs(title = paste("Price Index -", index_label),
         x = "Date", y = "Price index") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(size = 30, face = "bold"))
  
  # Put the charts together
  p1 + p2
}

plot_ita_ger <- plot_pair("ITA", "FTSEMIB") / plot_pair("GER", "DAX")
plot_fra_spa <- plot_pair("FRA", "CAC") / plot_pair("SPA", "IBEX")

# Save both charts separately
ggsave(file.path(base_path, "logreturns_priceindex_ITA_GER.pdf"),
       plot_ita_ger, width = 12, height = 10)

ggsave(file.path(base_path, "logreturns_priceindex_FRA_SPA.pdf"),
       plot_fra_spa, width = 12, height = 10)


# Histograms
index_labels <- c(
  "GER" = "DAX",
  "FRA" = "CAC",
  "ITA" = "FTSEMIB",
  "SPA" = "IBEX"
)

plot_histogram_with_normal <- function(df, index_name, color, xlim = c(-0.15, 0.15), ylim = c(0, 40)) {
  data_subset <- filter(df, Index == index_name)
  mu <- mean(data_subset$Log_Return)
  sigma <- sd(data_subset$Log_Return)
  
  ggplot(data_subset, aes(x = Log_Return)) +
    geom_histogram(aes(y = ..density..), bins = 50, fill = color, alpha = 0.6, boundary = 0) +
    stat_function(
      fun = dnorm,
      args = list(mean = mu, sd = sigma),
      color = "black", size = 1, linetype = "solid"
    ) +
    geom_vline(xintercept = mu, color = "black", linetype = "dotted", linewidth = 0.8) +
    labs(
      title = index_labels[index_name],  # Titolo sintetico desiderato
      x = "Log-Returns", y = "Density"
    ) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      panel.grid.minor = element_blank()
    )
}

p1 <- plot_histogram_with_normal(X_long, "ITA", index_colors["ITA"])  # FTSEMIB
p2 <- plot_histogram_with_normal(X_long, "GER", index_colors["GER"])  # DAX
p3 <- plot_histogram_with_normal(X_long, "FRA", index_colors["FRA"])  # CAC
p4 <- plot_histogram_with_normal(X_long, "SPA", index_colors["SPA"])  # IBEX

final_hist <- (p1 | p3) / (p2 | p4)

# Save pdf
ggsave(
  "Data/Chart preliminar analysis/histograms_logreturns_with_normal.pdf",
  final_hist,
  width = 12, height = 10
)



# QQ PLOT
X_long <- X %>%
  dplyr::select(ITA, GER, FRA, SPA) %>%
  pivot_longer(cols = everything(), names_to = "Index", values_to = "Value") %>%
  mutate(
    Index = recode(Index,
                   ITA = "FTSEMIB",
                   GER = "DAX",
                   FRA = "CAC",
                   SPA = "IBEX"),
    Index = factor(Index, levels = c("FTSEMIB", "CAC", "DAX", "IBEX"))  # ordine voluto
  ) %>%
  group_by(Index) %>%
  mutate(Value = scale(Value)) %>%
  ungroup()

qq_plot <- ggplot(X_long, aes(sample = Value)) +
  stat_qq(distribution = qnorm, size = 1, color = "blue") +
  stat_qq_line(distribution = qnorm, color = "red", linetype = "solid", linewidth = 0.8) +
  facet_wrap(~ Index, scales = "free", ncol = 2) +
  labs(
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(size = 30, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    plot.title = element_blank()
  )

ggsave(
  "Data/Chart preliminar analysis/qqplot_logreturns.pdf",
  qq_plot,
  width = 10, height = 8
)


# other Normality Tests

jarque.bera.test(X$ITA)
jarque.bera.test(X$GER)
jarque.bera.test(X$FRA)
jarque.bera.test(X$SPA)

shapiro.test(X$ITA)
shapiro.test(X$GER)
shapiro.test(X$FRA)
shapiro.test(X$SPA)

ad.test(X$ITA)
ad.test(X$GER)
ad.test(X$FRA)
ad.test(X$SPA)




# --------------------------------------------- Autocorrelation check ------------------------------------------------------------


# Funtion to improve basic ACF
plot_acf <- function(series, title) {
  ggAcf(series, lag.max = 30) +
    labs(title = title) +
    ylab(expression(Acf(r[t]))) +  # label y axis
    xlab("Lag") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 11)
    )
}

# Create the 4 plots
p1 <- plot_acf(X$ITA, "FTSEMIB")
p2 <- plot_acf(X$FRA, "CAC")
p3 <- plot_acf(X$GER, "DAX")
p4 <- plot_acf(X$SPA, "IBEX")

# Pair the 4 plots
acf_grid <- (p1 | p2) / (p3 | p4)

# Save in pdf
ggsave(
  "Data/Chart preliminar analysis/acf_logreturns.pdf",
  acf_grid,
  width = 12, height = 10
)


Box.test(X$ITA, type = "Ljung-Box", lag = 5)
Box.test(X$GER, type = "Ljung-Box", lag = 5)
Box.test(X$FRA, type = "Ljung-Box", lag = 5)
Box.test(X$SPA, type = "Ljung-Box", lag = 5)

Box.test(X$ITA, type = "Ljung-Box", lag = 10)
Box.test(X$GER, type = "Ljung-Box", lag = 10)
Box.test(X$FRA, type = "Ljung-Box", lag = 10)
Box.test(X$SPA, type = "Ljung-Box", lag = 10)

Box.test(X$ITA, type = "Ljung-Box", lag = 20)
Box.test(X$GER, type = "Ljung-Box", lag = 20)
Box.test(X$FRA, type = "Ljung-Box", lag = 20)
Box.test(X$SPA, type = "Ljung-Box", lag = 20)



# --------------------------- Volatility Clustering phenomenon ---------------------------------------

# ACF di |r_t|   (TAYLOR EFFECT)
plot_acf_abs <- function(series, title) {
  ggAcf(abs(series), lag.max = 30) +
    labs(title = title) +
    ylab(expression(Acf(abs(r[t])))) +
    xlab("Lag") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 11)
    )
}

p1 <- plot_acf_abs(X$ITA, "FTSEMIB")
p2 <- plot_acf_abs(X$FRA, "CAC")
p3 <- plot_acf_abs(X$GER, "DAX")
p4 <- plot_acf_abs(X$SPA, "IBEX")

acf_abs_grid <- (p1 | p2) / (p3 | p4)

ggsave(
  "Data/Chart preliminar analysis/acf_abs_logreturns.pdf",
  acf_abs_grid,
  width = 12, height = 10
)

#----

# Function for ACF r_t^2
plot_acf_squared <- function(series, title) {
  ggAcf(series^2, lag.max = 30) +
    labs(title = title) +
    ylab(expression(Acf(r[t]^2))) +
    xlab("Lag") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 11)
    )
}

p1 <- plot_acf_squared(X$ITA, "FTSEMIB")
p2 <- plot_acf_squared(X$FRA, "CAC")
p3 <- plot_acf_squared(X$GER, "DAX")
p4 <- plot_acf_squared(X$SPA, "IBEX")

acf_squared_grid <- (p1 | p2) / (p3 | p4)

ggsave(
  "Data/Chart preliminar analysis/acf_squared_logreturns.pdf",
  acf_squared_grid,
  width = 12, height = 10
)



Box.test(X$ITA^2, lag = 5, type = "Ljung-Box")
Box.test(X$GER^2, lag = 5, type = "Ljung-Box")
Box.test(X$FRA^2, lag = 5, type = "Ljung-Box")
Box.test(X$SPA^2, lag = 5, type = "Ljung-Box")

Box.test(X$ITA^2, lag = 10, type = "Ljung-Box")
Box.test(X$GER^2, lag = 10, type = "Ljung-Box")
Box.test(X$FRA^2, lag = 10, type = "Ljung-Box")
Box.test(X$SPA^2, lag = 10, type = "Ljung-Box")

Box.test(X$ITA^2, lag = 20, type = "Ljung-Box")
Box.test(X$GER^2, lag = 20, type = "Ljung-Box")
Box.test(X$FRA^2, lag = 20, type = "Ljung-Box")
Box.test(X$SPA^2, lag = 20, type = "Ljung-Box")





#---------------------------------------- MULTIVARIATE analysis ---------------------------------------------
# --------------------------------------- Correlation Matrix ------------------------------------------------
X_named <- X %>%
  dplyr::rename(
    FTSEMIB = ITA,
    DAX     = GER,
    CAC     = FRA,
    IBEX    = SPA
  ) %>%
  dplyr::select(FTSEMIB, DAX, CAC, IBEX)  # solo queste 4 colonne

# Pearson correlation matrix
cor_matrix_pearson <- cor(X_named, method = "pearson")

# Kendall correlation matrix
cor_matrix_kendall <- cor(X_named, method = "kendall")

# Plot Pearson
p_pearson <- ggcorrplot(
  cor_matrix_pearson,
  lab = TRUE,
  lab_size = 5,
  method = "square",
  type = "lower",
  colors = c("blue", "white", "red"),
  title = "",
  ggtheme = theme_minimal(base_size = 13)
) +
  theme(
    plot.title = element_blank(),
    axis.title = element_blank()
  )

# Plot Kendall
p_kendall <- ggcorrplot(
  cor_matrix_kendall,
  lab = TRUE,
  lab_size = 5,
  method = "square",
  type = "lower",
  colors = c("blue", "white", "red"),
  title = "",
  ggtheme = theme_minimal(base_size = 13)
) +
  theme(
    plot.title = element_blank(),
    axis.title = element_blank()
  )

ggsave("Data/Chart preliminar analysis/correlation_matrix_pearson_ggcorrplot.pdf",
       p_pearson, width = 6.5, height = 6)

ggsave("Data/Chart preliminar analysis/correlation_matrix_kendall_ggcorrplot.pdf",
       p_kendall, width = 6.5, height = 6)



# ----------------------- Lagged AutoCorrelation--------------------------

# labels
index_labels <- c("ITA" = "FTSEMIB", "GER" = "DAX", "FRA" = "CAC", "SPA" = "IBEX")
pairs <- combn(names(index_labels), 2, simplify = FALSE)

# Function CCF by ggAcf
plot_ccf_line <- function(series1, series2, label1, label2, squared = FALSE) {
  x <- X[[series1]]
  y <- X[[series2]]
  if (squared) {
    x <- x^2
    y <- y^2
  }
  
  ccf_obj <- ccf(x, y, lag.max = 20, plot = FALSE)
  df <- tibble(lag = ccf_obj$lag, acf = ccf_obj$acf)
  n <- length(x)
  conf_limit <- 2 / sqrt(n)
  
  ggplot(df, aes(x = lag, y = acf)) +
    geom_segment(aes(xend = lag, yend = 0), color = "grey30", linewidth = 0.6) +  # più sottile
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    geom_hline(yintercept = c(conf_limit, -conf_limit), linetype = "dashed", color = "blue", linewidth = 0.6) +
    labs(
      title = paste(label1, "vs", label2),
      x = "Lag", y = expression(Corr)
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 11)
    )
}

# Make the plots
plots_ccf  <- map(pairs, ~ plot_ccf_line(.x[1], .x[2], index_labels[.x[1]], index_labels[.x[2]], squared = FALSE))
plots_ccf2 <- map(pairs, ~ plot_ccf_line(.x[1], .x[2], index_labels[.x[1]], index_labels[.x[2]], squared = TRUE))

# put them in the same ppanel
ccf_plot  <- wrap_plots(plots_ccf,  ncol = 3)
ccf2_plot <- wrap_plots(plots_ccf2, ncol = 3)

ggsave("Data/Chart preliminar analysis/ccf_pairs_logreturns_lines.pdf",  ccf_plot,  width = 14, height = 10)
ggsave("Data/Chart preliminar analysis/ccf_pairs_squared_logreturns_lines.pdf", ccf2_plot, width = 14, height = 10)




# ------------------------ Tail dependence Empirical ------------------------------------


# labels
index_labels <- c("ITA" = "FTSEMIB", "GER" = "DAX", "FRA" = "CAC", "SPA" = "IBEX")

# Empirical Pseudo observation through ECDF
u <- map_dfc(X, ~ rank(.) / (length(.) + 1)) %>%
  dplyr::select(ITA, GER, FRA, SPA)

# Create the scatterplots and the blue and red lines
pseudo_pairs <- combn(names(u), 2, simplify = FALSE)
plot_list <- map(pseudo_pairs, function(pair) {
  df <- tibble(
    x = u[[pair[1]]],
    y = u[[pair[2]]]
  ) %>%
    mutate(
      lower_tail = (x < 0.05) & (y < 0.05),
      upper_tail = (x > 0.95) & (y > 0.95)
    )
  
  ggplot(df, aes(x = x, y = y)) +
    geom_point(data = filter(df, !lower_tail & !upper_tail), color = "gray50", alpha = 0.6, size = 1.2) +
    geom_point(data = filter(df, lower_tail), color = "red", size = 1.2) +
    geom_point(data = filter(df, upper_tail), color = "blue", size = 1.2) +
    geom_hline(yintercept = c(0.05, 0.95), linetype = "dashed", color = c("red", "blue")) +
    geom_vline(xintercept = c(0.05, 0.95), linetype = "dashed", color = c("red", "blue")) +
    coord_equal() +
    labs(
      x = index_labels[pair[1]],
      y = index_labels[pair[2]]
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_blank(),  # no title
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 10)
    )
})

pseudo_grid <- wrap_plots(plot_list, nrow = 2, ncol = 3)

ggsave(
  "Data/Chart preliminar analysis/pseudo_observations_pairs.pdf",
  pseudo_grid,
  width = 12, height = 12
)
