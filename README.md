# Comparison Of GARCH-Copula And An Optimal Transport-Based Method

This repository contains the code, data, and materials related to my master's thesis in Finance, completed as part of the Double Degree program between the University of Pavia (Italy) and HEC Liège (Belgium).

## Abstract
[Insert here the full abstract of your thesis.]

## Data Source
- **Market data**: Daily closing prices for FTSEMIB, DAX40, CAC40, and IBEX35 from January 5, 2015 to April 25, 2025.  
- **Provider**: Investing.com.  
- **Frequency**: Daily observations, adjusted for corporate actions (splits, dividends).  

## Installation & Requirements

This project uses both **R** and **MATLAB**.  
R is used for the GARCH-Copula modelling, backtesting, and results analysis.  
MATLAB is used for the Optimal Transport (OT) approach with the 4D Sinkhorn algorithm.

---

### R Environment
- **Version**: R (≥ 4.4.2)  
- **Required packages**:  
  `dplyr`, `ggplot2`, `xts`, `tidyr`, `moments`, `tseries`, `forecast`, `MASS`,  
  `copula`, `rugarch`, `nortest`, `fitdistrplus`, `stats`, `FinTS`, `GGally`,  
  `gridExtra`, `profvis`, `parallel`, `patchwork`, `zoo`, `tidyverse`,  
  `ggcorrplot`, `grid`, `knitr`, `kableExtra`

To install all R dependencies at once:
```r
install.packages(c(
  "dplyr","ggplot2","xts","tidyr","moments","tseries","forecast","MASS",
  "copula","rugarch","nortest","fitdistrplus","stats","FinTS","GGally",
  "gridExtra","profvis","parallel","patchwork","zoo","tidyverse","ggcorrplot",
  "grid","knitr","kableExtra"
))

### MATLAB Environment
- **Version**: MATLAB R2018b or later (recommended ≥ R2022a)  
- **Required Toolbox**: Statistics and Machine Learning Toolbox  
  (needed for `randsample`, `corr('type','Kendall')`, `skewness`, and `kurtosis`)




## Project Structure
The repository is organized as follows:

