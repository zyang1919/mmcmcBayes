# mmcmcBayes

**mmcmcBayes** is an R package that implements **Multistage Markov Chain Monte Carlo (MMCMC) method** for detecting **Differentially Methylated Regions (DMRs)** 
in **DNA methylation** data. It uses Bayesian inference with **Alpha-Skew Generalized Normal (ASGN) model** and supports **Bayes Factor (BF)** or **Anderson-Darling (AD) test** for region selection.

---

## Features
- Implements **Multistage MMCMC** to refine DMR detection.
- Supports **two statistical testing approaches**:
  - **Bayes Factor (BF)** – Bayesian model comparison for DMR detection.
  - **Anderson-Darling Test (AD)** – Nonparametric frequentist test for distributional differences.
- Provides the location of detected DMRs using the CpG coordinates as **Start_CpG** and **End_CpG**.
- Provides **CpG_Count**, showing how many CpG sites detected in the region.
- Provides **Decision_Value** from the test methods. 
- Supports **customizable MCMC settings and priors**.



## Function Arguments

### `mmcmcBayes()` ###
- **`cancer_data`**, **`normal_data`**: Data frames containing **M-values** for cancer and normal groups.
- **`test`**: Statistical test for DMR detection (`"BF"` for Bayes Factor, `"AD"` for Anderson-Darling Test).
- **`max_stages`**: Maximum number of iterative refinement stages (default: `3`).
- **`num_splits`**: Number of sub-segments per stage (default: `10`).
- **`bf_thresholds`**: Bayes Factor thresholds per stage. If `NULL`, then use the default value (`list(stage1 = 10, stage2 = 15, stage3 = 20)`).
- **`pvalue`**: AD test p-value thresholds. If `NULL`, then use the default value (`list(stage1 = 1e-4, stage2 = 1e-6, stage3 = 1e-8)`).
- **`mcmc`**: List of MCMC settings (`list(nburn = 5000, niter = 10000, thin = 1)`).
- **`priors`**: List of prior values (default: `NULL`).

### `compare_dmrs()` ###
- **`rst1`**: Detected DMRs from one method.
- **`rst2`**: Detected DMRs from another method. 
---

## Installation
To install `mmcmcBayes`, **clone the repository** and **load the functions**:

```r
# Clone the repository
system("git clone https://github.com/zyang1919/mmcmcBayes.git")

# Load the functions
source("your_path/mmcmcBayes/mmcmcBayes.R")
```

## Examples
```r
## read data
wk <- "Your path for the data/Data"
data1 <- get(load(paste(wk,"chr6_methylation_M_values_cancer_retrospective.RData",sep='/')))
data2 <- get(load(paste(wk,"chr6_methylation_M_values_normal_retrospective.RData",sep='/')))

## This is an example of using Bayes Factor
result1 <- mmcmcBayes(data1, data2, stage = 1,max_stages = 3,num_splits = 10,
                   test = "BF",
                   bf_thresholds = list(stage1 = 10.6, stage2 = 10.7, stage3 = 10.8))

## This is an example of us Anderson-Darling Test
result2 <- mmcmcBayes(data1, data2, stage = 1,max_stages = 3,num_splits = 10,
                   test = "AD",
                   pvalue = list(stage1 = 10^(-6), stage2 = 10^(-6), stage3 = 10^(-6)))

## Comparing the results from two test
overlap <- compare_dmrs(result1,result2)

```


## Model Output
```r
result1[1:8,]

   Chromosome  Start_CpG    End_CpG CpG_Count Decision_Value
1           6 cg00000721 cg00201275       365       10.88718
2           6 cg00201779 cg00446211       364       10.93774
3           6 cg00669964 cg00944666       365       11.04935
4           6 cg00944712 cg01178040       364       10.80733
5           6 cg01180523 cg01414358       364       10.88814
6           6 cg01414663 cg01664325       365       10.90684
7           6 cg01664382 cg01916115       364       10.87119
8           6 cg01916632 cg02151997       364       10.82806

result2[1:9,]

   Chromosome  Start_CpG    End_CpG CpG_Count Decision_Value
1           6 cg00000721 cg00201275       365     2.8596e-09
2           6 cg00201779 cg00446211       364     2.9090e-09
3           6 cg00944712 cg01178040       364     5.2034e-07
4           6 cg01180523 cg01414358       364     2.4175e-08
5           6 cg01414663 cg01664325       365     1.6910e-07
6           6 cg01664382 cg01916115       364     3.8775e-10
7           6 cg01916632 cg02151997       364     6.4501e-08
8           6 cg02152351 cg02407730       365     1.0132e-07
9           6 cg02407762 cg02689448       365     3.3667e-12

overlap[1:9,]

   Chromosome Start_CpG_Method1 End_CpG_Method1 Start_CpG_Method2 End_CpG_Method2 Overlap_Percentage
1           6        cg00000721      cg00201275        cg00000721      cg00201275                100
2           6        cg00201779      cg00446211        cg00201779      cg00446211                100
3           6        cg00944712      cg01178040        cg00944712      cg01178040                100
4           6        cg01180523      cg01414358        cg01180523      cg01414358                100
5           6        cg01414663      cg01664325        cg01414663      cg01664325                100
6           6        cg01664382      cg01916115        cg01664382      cg01916115                100
7           6        cg01916632      cg02151997        cg01916632      cg02151997                100
8           6        cg02152351      cg02407730        cg02152351      cg02407730                100
9           6        cg02407762      cg02689448        cg02407762      cg02689448                100
```

### **Authors**
This package was developed by:

- **Zhexuan Yang**, Ph.D., Northern Illinois University   
- **Duchwan Ryu**, Ph.D., Northern Illinois University  
- **Feng Luan**, Ph.D. Candidate, Northern Illinois University

