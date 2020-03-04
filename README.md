# single_index

This repository gives R scripts for simulations with the ordinary profile Least Squares Estimator (LSE), the Simple Score Estimator (SSE), the Effective Score Estimator (ESE). It also gives simulations for:
1. a spline estimator (Arun K. Kuchibhotla and Rohit K.Patra. Efficient estimation in single index models through smoothing splines. Bernoulli, 26(2):1587–1618, 2020, https://doi.org/10.3150/19-BEJ1183). This is implemented in the R package simest, but we use our own implementation.
2. the EDR estimate (Marian Hristache, Anatoli Juditsky, and Vladimir Spokoiny. Direct estimation of the index coefficient in a single-index model. Ann. Statist., 29(3):595–623, 2001. ISSN 0090-5364. https://doi.org/10.1214/aos/1009210681),
3. the MAVE method (from the R package, theory is in Yingcun Xia. Asymptotic distributions for two estimators of the single-index model. Econometric Theory, 22(6):1112–1137, 2006. ISSN 0266-4666, https://doi.org/10.1017/S0266466606060531)
4. the EFM algorithm (Xia Cui, Wolfgang Karl Härdle, and Lixing Zhu. The EFM approach for single-index models. Ann. Statist., 39(3):1658–1688,  2011, https://doi.org/10.1214/10-AOS871).

The methods are discussed in "Profile least squares estimators in the monotone single index model (2020)" https://arxiv.org/abs/2001.05454 of Fadoua Balabdaoui and Piet Groeneboom.

The R script simulation_table1.R in the directory "Comparisons" runs simulations for the first model, discussed in "Profile least squares estimators in the monotone single index model (2020)" and the R script simulation_table2.R runs simulations for the second model. One can adjust the sample size n and the number of replications NumIt in these R files.
For running the scripts, one needs to be able to run the R package Rcpp. The values of the simulations are put into the growing files LSE.txt, etc. If one wants to have new files for a new run, one has to remove the previous files of this type. The means and covariances of the runs can be computed by the script variances.R.

