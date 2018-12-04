# single_index

This repository gives Rcpp scripts for simulations with the profile Least Squares Estimator (LSE), the Simple Score Estimator (SSE), the Effective Score Estimator (ESE), the Penalized Score Estimator (PLSE),
a linear estimate and Han's Maximum Rank Correlation Estimate, using the methods of section 4.2 and 5 of "Score estimation in the monotone single index model" by Fadoua Balabdaoui, Piet Groeneboom and Kim Hendrickx, https://doi.org/10.1111/sjos.12361, Scandinavian Journal of Statistics and also the manuscript ``The Lagrange approach in the monotone single index model'', http://dutiosb.twi.tudelft.nl/%7Epietg/Single_index_Lagrange.pdf of Piet Groeneboom.

The R script simulation_table2.R (for Table 2 in the paper) in the directory "Comparisons"
compares the simple score estimator (SSE), the efficient score estimator (ESE) and the least squares estimator (LSE) with the algorithm of the R package EDR and the penalized least squares estimator (PLSE), using smoothing splines.

I could not remove the text that is produced during the runs for the EDR package.
In the diectory ``Comparisons'' the PLSE is computed using the R package "simest", with the slight modification "NewSimEst.R", kindly provided to me by the authors of the package. At the end of the run, the least squares estimate of the link function, corresponding to the ESE algorithm, is shown for the last sample (blue step function), together with the underlying link function (red dashed function) and the points of the sample, plotted in the scale provided by the estimate of the regression parameter alpha.

In the directory ``Single_index_Lagrange'' the PLSE is computed by the Lagrange method, and here the PLSE, computed using the R package "simest", is called PLSE2. This directory is the companion of the manuscript http://dutiosb.twi.tudelft.nl/%7Epietg/Single_index_Lagrange.pdf.

For running the scripts, one needs to be able to run the package Rcpp. One also needs to
copy the other files of the directory, since most of them are used in the script.

