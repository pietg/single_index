# single_index

This repository gives R scripts for simulations with the ordinary profile Least Squares Estimator (LSE), the Simple Score Estimator (SSE), the Effective Score Estimator (ESE) and a spline estimator, using the methods explained in "Profile least squares estimators in the monotone single index model (2020)" https://arxiv.org/abs/2001.05454 of Fadoua Balabdaoui and Piet Groeneboom.

The R scripts in the directory "Comparisons" compare the ordinary profile least squares estimator (LSE), the simple score estimator (SSE) and the efficient score estimator (ESE) with the algorithm of the R package EDR and a spline estimator. We use our own implementation of the spline estimator, see https://arxiv.org/abs/2001.05454.
The augmented Lagrange method is used to avoid reparametrization, in combination with the Hooke-Jeeves pattern search algorithm, see again https://arxiv.org/abs/2001.05454.

For running the scripts, one needs to be able to run the package Rcpp. One also needs to
copy the other files of the directory, since most of them are used in the script. We give as examples the R scripts for 1000 replications of sample size 100 for the two models, used in the simulations in "Profile least squares estimators in the monotone single index model (2020)" https://arxiv.org/abs/2001.05454. Simulation1.R is for the first model, with psi(x)=x^3 and the error independent of the covariates. Simulation2.R is for the logistic-binomial model, with an error that depends on the covariates. The values of the simulation are put into the growing files LSE.txt, etc. If one wants to have new files for a new run, one has to remove the previous files of this type. The means and covariances of the runs can be computed by the script table.R.

