# SI-FARM: Confounder adjustment in single index function-on-scalar regression model

# Introduction
## Abstract
The function-on-scalar regression model serves as a potent tool for elucidating the connection between functional responses and covariates
 of interest. Despite its widespread utilization in numerous extensive neu
roimaging investigations, prevailing methods often fall short in account
ing for the intricate nonlinear relationships and the enigmatic confound
ing factors stemming from imaging heterogeneity. This heterogeneity may
 originate from a myriad of sources, such as variations in study environ
ments, populations, designs, protocols, and concealed variables. To address
 this challenge, this paper develops a single index function-on-scalar regres
sion model to investigate the nonlinear associations between functional re
sponses and covariates of interest while making adjustments for concealed
 confounding factors arising from potential imaging heterogeneity. Both es
timation and inference procedures are established for unknown parameters
 within our proposed model. In addition, the asymptotic properties of esti
mated functions and detected confounding factors are also systematically
 investigated. The finite-sample performance of our proposed method is as
sessed by using both Monte Carlo simulations and a real data example on
 the diffusion tensor images from the Alzheimer’s Disease Neuroimaging Ini
tiative study.

## Paper
[Ding, S., Zhou, X., Lin, J., Liu, R., & Huang, C. (2024). Confounder adjustment in single index function-on-scalar regression model. Electronic Journal of Statistics, 18(2), 5679-5714.](https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-18/issue-2/Confounder-adjustment-in-single-index-function-on-scalar-regression-model/10.1214/24-EJS2333.full)

 <br>

# Folder Structure
- **[./utilities/](./utilities/):** 
	contains all the user-defined functions
- **[./SIVC_2016_code](./SIVC_2016_code/):**
    the refernece code for comparison with [Luo, X., Zhu, L., & Zhu, H. (2016).](https://onlinelibrary.wiley.com/doi/full/10.1111/biom.12526)

<br>


# Code Description
## [simu_main.m](simu_main.m)
> **the script for calculating the estimation /prediction errors from simulated datasets measured by the mean and standard deviation (std) of ISE /SE.**
### (1) Settings
> **For different simulation settings, please modify the corresponding parameters.**
- n: sample size, choose from [100,200]
- p: number of covariates, p = 3
- nv: number of grids for functions, choose from [100,200]
- nsimu: number of simulated datasets, nsimu = 100
- ul: the level of heterogeneity, choose from [0,1,2,3], corresponding to np, weak, moderate, and high heterogeneity respectively
- index_g: the choice of index function, choose from [1,2], corresponding to (1) g(t)=sin(2t)+2cos(2+t) and (2) g(t)=exp(−t) respectively

### (2) An example 
> **Settings: n=100, p=3, nv=200, nsimu=100, ul=3, index_g=1.**
**Run the command through MATLAB command prompt**
```matlab
simu_main(100, 3, 200, 100, 3, 1);
```

<br>

## [factor_analysis.R](factor_analysis.R)
> **the script for factor analysis to determine which method to use, then obtian q_est by the chosen method.**
> ### (1) Settings
> **For different simulation settings, please modify the corresponding parameters.**
- n: sample size, choose from [100,200]
- p: number of covariates, p = 3
- nv: number of grids for functions, choose from [100,200]
- nsimu: number of simulated datasets, nsimu = 100
- ul: the level of heterogeneity, choose from [0,1,2,3], corresponding to np, weak, moderate, and high heterogeneity respectively
- index_g: the choice of index function, choose from [1,2], corresponding to (1) g(t)=sin(2t)+2cos(2+t) and (2) g(t)=exp(−t) respectively

### (2) An example 
> **Settings: n=100, p=3, nv=200, nsimu=100, ul=3, index_g=1.**
**Run the R script "factor_analysis.R" using the command line:**
```r
Rscript factor_analysis.R 100 3 200 100 3 1
```

## [simu_coverage_prob.m](simu_coverage_prob.m)
> **the script for calculating the coverage probabilities for SCB of functional coefficients $\beta(s)$ & $g(x\beta(s))$**

<br>
<br>


# References
- [Luo, X., Zhu, L., & Zhu, H. (2016). Single‐index varying coefficient model for functional responses. Biometrics, 72(4), 1275-1284.](https://onlinelibrary.wiley.com/doi/full/10.1111/biom.12526)
