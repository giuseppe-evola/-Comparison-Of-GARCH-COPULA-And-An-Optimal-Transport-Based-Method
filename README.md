# Comparison Of GARCH-Copula And An Optimal Transport-Based Method

This repository contains the code, data, and materials related to my master's thesis in Finance. I am pursuing my degree within the Double Degree program between the University of Pavia (Italy) and HEC Liège (Belgium).

## Abstract
[Insert here the full abstract of your thesis.]

## Data Source
- **Market data**: Daily closing prices for FTSEMIB, DAX40, CAC40, and IBEX35 from January 5, 2015 to April 25, 2025.  
- **Provider**: Investing.com.  
- **Frequency**: Daily observations, adjusted for corporate actions (splits, dividends).  

## Installation & Requirements

This project uses both **R** and **MATLAB**.  
R is used for the GARCH and Copula models definition and GARCH-Copula and GARCH-OT modelling for VaR estimation.  
MATLAB is used for the application of the Optimal Transport-based model.

### R Environment  
**Version**: R (≥ 4.4.2)  
**Required packages**:  
`dplyr`, `ggplot2`, `xts`, `tidyr`, `moments`, `tseries`, `forecast`, `MASS`,  
`copula`, `rugarch`, `nortest`, `fitdistrplus`, `stats`, `FinTS`, `GGally`,  
`gridExtra`, `profvis`, `parallel`, `patchwork`, `zoo`, `tidyverse`,  
`ggcorrplot`, `grid`, `knitr`, `kableExtra`  

**To install all R dependencies at once:**  
 ```r
 install.packages(c(
   "dplyr","ggplot2","xts","tidyr","moments","tseries","forecast","MASS",
   "copula","rugarch","nortest","fitdistrplus","stats","FinTS","GGally",
   "gridExtra","profvis","parallel","patchwork","zoo","tidyverse","ggcorrplot",
   "grid","knitr","kableExtra"
 ))
 ```

### MATLAB Environment  
**Version**: MATLAB R2018b or later (recommended ≥ R2022a)  
**Required Toolbox**: Statistics and Machine Learning Toolbox  
(needed for `randsample`, `corr('type','Kendall')`, `skewness`, and `kurtosis`)

---

## **Project Structure**

The project is organized into four main sections:

### **1. Preliminary Analysis**
- **`1_preliminary_analysis.R`**  
  Splits the dataset into training and test sets, computes descriptive statistics, and highlights the stylized facts observed in the analyzed time series.

### **2. Model Selection and Functions**
- **`Functions.R`**  
  Contains all custom functions used across the project, including model selection procedures and transformation tools.
- **`2_best_garch_copula.R`**  
  - **Part 1**: Identifies the optimal ARMA-GARCH model for each stock index.  
  - **Part 2**: Performs a preliminary analysis of the best copula based on standardized residuals.

### **3. GARCH-Copula Approach**
- **`3.1_GARCH_COPULA_DYNAMIC.R`**  
  Implements the rolling-window GARCH-Copula algorithm with dynamic refitting for VaR estimation. Parallelization is used to improve computational efficiency. Includes visualizations of the time evolution of both marginal and copula parameters.
- **`3.1_GARCH_COPULA_STATIC.R`**  
  Implements the static GARCH-Copula approach, where copula models are estimated only on the training set and applied to all rolling windows.

### **4. GARCH-OT (MATLAB + R)**

**In the `OT Algorithm definitive` directory:**
- **`sinkhorn_mot_4.m`** (MATLAB)  
  Core implementation of the regularized multimarginal Optimal Transport algorithm using the Sinkhorn method.
- **`algo_parameters_test.m`** (MATLAB)  
  Executes the Sinkhorn algorithm across 16 different hyperparameter combinations (ε, n) and selects the top 3 using a custom selection procedure.
- **`ToBeUsed.m`** (MATLAB)  
  Main script that generates the OT-based simulations (three 100000x4 datasets) from the top three selected specifications (ε, n).

**In the main directory:**
- **`4_GARCH_OT.R`**  
  Loads the top three OT-based simulations and applies the GARCH-OT algorithm to estimate the VaR. **Note**: these simulations are of standardized residuals, unlike copula-based approaches which use pseudo-observations.

### **Other**
- **`other.R`**  
  Contains auxiliary scripts used throughout the thesis.

---

## **How to Reproduce the Analysis**

A `.RData` file containing all variables and results is included in the repository, allowing you to quickly reproduce the thesis outcomes without re-running the entire pipeline.

If you want to execute the full analysis from scratch, follow these steps:

**0.** Download the entire repository and set your working directory.

**1.** Run **`1_preliminary_analysis.R`**.

**2.** Run **`2_best_garch_copula.R`**.  
This script will export the standardized residuals obtained using the selected GARCH models on the training set. The output must be manually moved into the `OT Algorithm definitive` directory (a placeholder file is already present there). This choice is intentional to maintain full control over the file, although automatic saving could be implemented by modifying the output path.

### **For GARCH-Copula approaches and static vs dynamic comparison:**
**a.** Run **`3.1_GARCH_COPULA_DYNAMIC.R`**  
**b.** Run **`3.1_GARCH_COPULA_STATIC.R`**  

### **For the GARCH-OT approach:**
**a.** In the `OT Algorithm definitive` directory, run **`ToBeUsed.m`** (MATLAB).  
This will produce three CSV files containing OT-based simulations for the top three specifications:  
- `sinkhorn_simulations_TOP1.csv`  
- `sinkhorn_simulations_TOP2.csv`  
- `sinkhorn_simulations_TOP3.csv`  

**b.** Run **`4_GARCH_OT.R`** in the main directory.















