# Comparison Of GARCH-Copula And An Optimal Transport-Based Method

This repository contains the code, data, and materials related to my master's thesis in Finance. I am pursuing my degree within the Double Degree program between the University of Pavia (Italy) and HEC Liège (Belgium).

## Abstract
Accurately estimating the Value at Risk (VaR) constitutes a central challenge both in
the academic literature and in financial practice. When the problem is analyzed from
a multivariate perspective, the main difficulty lies in appropriately capturing the depen-
dence structure among different assets. This thesis investigates, following an agnostic
approach, two alternative methods to address this problem in the context of portfolio mul-
tivariate VaR estimation, both based on GARCH-type models to describe the univariate
dynamics of the asset returns. The first approach, of parametric nature, relies on cop-
ula functions (GARCH-Copula approach); the second constitutes the innovative aspect of
this work, exploiting Optimal Transport theory to model dependence in a non-parametric
way, thus leading to the GARCH-OT approach. The performance of GARCH-Copula
models is assessed using the six most common static copulas (Gaussian, t, Frank, Gum-
bel, Clayton, and Joe), and compared with those of the three optimal specifications of the
GARCH-OT model, calibrated through a novel procedure developed in this study. The
analysis employs a rolling window framework, estimating VaR at the 95% and 99% con-
fidence levels for an equally weighted portfolio composed of four European stock indices,
and the performances are evaluated using the Kupiec test, the Christoffersen test, and the
Lopez statistic. The results highlight the effectiveness of both approaches for VaR esti-
mation. In the GARCH-Copula context, models based on elliptical copulas outperform
Archimedean ones, with the exception of the Clayton copula, whose performance is com-
parable to that of elliptical models. The GARCH-OT represents a promising alternative,
yielding comparable results at the 99% level and slightly underestimating risk at the 95%
level relative to the GARCH-Copula, suggesting significant potential, although requiring
higher computational cost and accurate calibration.

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

**1.** Run **`1_preliminary_analysis.R`** and **`Functions.R`**

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















