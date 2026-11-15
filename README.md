[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3936045.svg)](https://doi.org/10.5281/zenodo.3936045)

# longCombat: Longitudinal ComBat R Package

Longitudinal ComBat uses an empirical Bayes method to harmonize means and variances of the residuals across batches in a linear mixed effects model framework. Detailed methods are described in the manuscript: 

Beer JC, Tustison NJ, Cook PA, Davatzikos C, Sheline YI, Shinohara RT, Linn KA. (2020) Longitudinal ComBat: A method for harmonizing longitudinal multi-scanner imaging data. NeuroImage. In press. https://doi.org/10.1016/j.neuroimage.2020.117129.

Install package with: 
```{r, include=FALSE}
install.packages("devtools")
devtools::install_github("jcbeer/longCombat")
```

Note: longCombat currently will not run if tidyverse suite is loaded. This may be fixed in the future. For now, please run longCombat before loading tidyverse.

[Contact Joanne Beer](mailto:joannecbeer@gmail.com?subject=[GitHub]%20longCombat) with any questions, comments, or suggestions. Feedback is appreciated. 

# Applying longComBat to a test set.

The file LongComBat_apply_test.R contains an examplpe of how to apply the parameters from longgitudinal ComBat to a test set.
The steps are the following:

1. Apply a modified version of longComBat to the training set. The only modification is that the file contains a container to extract the fixed effects and variance components (between subjects variance and within subject variance).

2. Remove the overall intercept alpha, the covariates estimates and the site estimates from the test set. What remains is the residuals per scan. Sum those residuals per subject. Then, instead of averaging per subject, divide by a shrunken average eta_j  <- sum_e / (n_j + sigma^2 / tau2) (adding a factor taking into account the ratio of between and within subject noise)

3. Standardize, remove site and subject factors and propagate back.
