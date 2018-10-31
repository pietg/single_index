# single_index

This repository gives Rcpp scripts for simulations with the simple score estimator and a
link-free estimate, using the methods of section 4.2 and 5 of "Score estimation in the
monotone single index model" by Fadoua Balabdaoui, Piet Groeneboom and Kim Hendrickx,
http://dutiosb.twi.tudelft.nl/%7Epietg/BGH_final.pdf, to appear in the Scandinavian Journal of Statistics.

The R script simulation_table2.R (for Table 2 in the paper) in the directory "Comparisons"
compares the simple score estimator (SSE), the efficient score estimator (ESE) and the least squares estimator (LSE) with the algorithm of the R package EDR and the penalized least squares estimator (PLSE), using smoothing splines.

I could not remove the text that is produced during the runs for the EDR package.
The PLSE is computed using the R package "simest", with the slight modification "NewSimEst.R", kindly provided to me by the authors of the package. At the end of the run, the least squares estimate of the link function, corresponding to the ESE algorithm, is shown for the last sample (blue step function), together with the underlying link function (red dashed function) and the points of the sample, plotted in the scale provided by the estimate of the regression parameter alpha.

For running the scripts, one needs to be able to run the package Rcpp. One also needs to
copy the other files of the directory, since most of them are used in the script.
