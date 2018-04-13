Feature Selection Hidden Markov Model
===========

R-package that performs a feature selection strategy that selects the best variables of the observation matrix of an experiment with one or multiple conditions by using a hidden Markov model. This strategy can summarize the samples observations or receive a priori information to get a customized model. The multiple condition feature selection strategy depends on the Control/Baseline condition to have a more accurate filtering.

Quick setup
-----------
- Install R on your machine
- Install dependencies from CRAN repository:
  - Rcpp, RcppArmadillo, RcppHMM
- Install the source package FSHMM_1.0.3.tar.gz
  - From R console:
  ```R
  setwd("PATH") # point this to the source file location
  install.packages("./FSHMM_1.0.2.tar.gz", repos = NULL, type = "source")
  ```
- Load library
```R
library(FSHMM)
```

Main functions
-----------
The functions used for feature selection are:
-featureSelectionHMM for case-control or multiple condition feature selection
-featureSelectionOneConditionHMM for unique condition feature selection

The functions used for plotting the top best ranked features
-plotMultivariate for studies where all the replicates are used for model parameter estimation
-plotUnivariate for studies where all the replicates are summarized for model parameter estimation


