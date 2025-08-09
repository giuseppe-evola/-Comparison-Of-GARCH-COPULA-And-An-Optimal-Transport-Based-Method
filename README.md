# Comparison Of GARCH-Copula And An Optimal Transport-Based Method

This repository contains the code, data, and materials related to my master's thesis in Finance, completed as part of the Double Degree program between the University of Pavia (Italy) and HEC Liège (Belgium).

## Abstract
[Insert here the full abstract of your thesis.]

## Data Source
- **Market data**: Daily closing prices for FTSEMIB, DAX40, CAC40, and IBEX35 from January 1, 2015 to April 26, 2025.  
- **Provider**: Investing.com.  
- **Frequency**: Daily observations, adjusted for corporate actions (splits, dividends).  

## Installation & Requirements
- **Programming language**: R (version ≥ 4.3.0)  
- **Main dependencies**:  
  - `rugarch`  
  - `copula`  
  - `parallel`  
  - `ggplot2`  
  - `dplyr`  

To install all dependencies at once:
```r
install.packages(c("rugarch", "copula", "parallel", "ggplot2", "dplyr"))
