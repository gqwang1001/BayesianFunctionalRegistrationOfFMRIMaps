#### This repo provides the code for the paper `Bayesian Functional Registration of fMRI Activation Maps`
Please read the instructions below for better understanding.

1. Perform the corresponding features-based prior estimation through the `priors_calculation_symmetric.m`. please adjust the code for your need.
2. Once you have the estimation of parameters for the prior distribution, you are ready to run the sampling in STAN through `sym_main.m` file.
3. Choose the best hyperparameters with PSIS-LOO through `Loo_comparison_sym.R`. 
4. Compute the credible regions in `CredibleRegions.R`.





