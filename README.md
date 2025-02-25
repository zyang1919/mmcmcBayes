# mmcmcBayes

**mmcmcBayes** is an R package that implements **Multistage Markov Chain Monte Carlo (MMCMC) method** for detecting **Differentially Methylated Regions (DMRs)** 
in **DNA methylation** data. It uses Bayesian inference with **Alpha-Skew Generalized Normal (ASGN) model** and supports **Bayes Factor (BF)** or **Anderson-Darling (AD) test** for region selection.

---

## Features
- Implements **Multistage MMCMC** to refine DMR detection.
- Supports **two statistical testing approaches**:
  - **Bayes Factor (BF)** – Bayesian model comparison for DMR detection.
  - **Anderson-Darling Test (AD)** – Non-parametric Frequentist test for distributional differences.
- Provides the location of detected DMRs using the CpG coordinates as **Start_CpG** and **End_CpG**.
- Provides **CpG_Count**, showing how many CpG sites detected in the region.
- Provides **Decision_Value** from the test methods. 
- Supports **customizable MCMC settings and priors**.

---

## Function Arguments

### **Main Function: `mmcmcBayes()`**
- **`cancer_data`**, **`normal_data`**: Data frames containing **M-values** for cancer and normal groups.
- **`test`**: Statistical test for DMR detection (`"BF"` for Bayes Factor, `"AD"` for Anderson-Darling Test).
- **`max_stages`**: Maximum number of iterative refinement stages (default: `3`).
- **`num_splits`**: Number of sub-segments per stage (default: `10`).
- **`bf_thresholds`**: Bayes Factor thresholds per stage. If `NULL`, then use the default value (`list(stage1 = 10, stage2 = 15, stage3 = 20)`).
- **`pvalue_thresholds`**: AD test p-value thresholds. If `NULL`, then use the default value (`list(stage1 = 1e-4, stage2 = 1e-6, stage3 = 1e-8)`).
- **`mcmc`**: List of MCMC settings (`list(nburn = 5000, niter = 10000, thin = 1)`).
- **`priors`**: List of prior values (default: `NULL`).

---

## Installation
To install `mmcmcBayes`, **clone the repository** and **load the functions**:

```r
# Clone the repository
system("git clone https://github.com/zyang1919/mmcmcBayes.git")

# Load the functions
source("your_path/mmcmcBayes/mmcmcBayes.R")
```

--

## Example
```r
wk <- "Your path for the data/Data"
data1 <- get(load(paste(wk,"chr6_methylation_M_values_cancer_retrospective.RData", 
                        sep='/')))
data2 <- get(load(paste(wk,"chr6_methylation_M_values_normal_retrospective.RData", 
                        sep='/')))

## This is an example of using Bayes Factor
result <- mmcmcBayes(data1, data2, stage = 1,max_stages = 3,num_splits = 10,
                   test = "BF",
                   bf_thresholds = list(stage1 = 10.6, stage2 = 10.7, stage3 = 10.8))
```



