# Comparison Of GARCH-Copula And An Optimal Transport-Based Method

This repository contains the code, data, and materials related to my master's thesis in Finance, completed as part of the Double Degree program between the University of Pavia (Italy) and HEC Liège (Belgium).

## Abstract
[Insert here the full abstract of your thesis.]

## Data Source
- **Market data**: Daily closing prices for FTSEMIB, DAX40, CAC40, and IBEX35 from January 5, 2015 to April 25, 2025.  
- **Provider**: Investing.com.  
- **Frequency**: Daily observations, adjusted for corporate actions (splits, dividends).  

# Installation & Requirements

## R Environment
- **R** (Version >= 4.4.2)

### Required R Packages
The following R packages are required for the analysis:

`dplyr`, `ggplot2`, `xts`, `tidyr`, `moments`, `tseries`, `forecast`, `MASS`, `copula`, `rugarch`, `nortest`, `fitdistrplus`, `stats`, `FinTS`, `GGally`, `gridExtra`, `profvis`, `parallel`, `patchwork`, `zoo`, `tidyverse`, `ggcorrplot`, `grid`, `knitr`, `kableExtra`

### Installation
To install all dependencies at once, run the following command in R:

```r
install.packages(c(
 "dplyr", "ggplot2", "xts", "tidyr", "moments", "tseries", "forecast", "MASS",
 "copula", "rugarch", "nortest", "fitdistrplus", "stats", "FinTS", "GGally",
 "gridExtra", "profvis", "parallel", "patchwork", "zoo", "tidyverse", "ggcorrplot",
 "grid", "knitr", "kableExtra"
))



## Project Structure
The repository is organized as follows:

