# single_index

This repository gives Rcpp scripts for simulations with the simple score estimator and a
link-free estimate, using the methods of section 4.2 and 5 of "Score estimation in the
monotone single index model" by Fadoua Balabdaoui, Piet Groeneboom and Kim Hendrickx,
https://arxiv.org/abs/1712.05593

The R scripts simulation_table1.R (for Table 1 in the paper) and simulation_table2.R
(for Table 2 in the paper) in the directory ``Comparisons'' compare the simple score
estimator (SSE), the efficient score estimator (ESE) and the least squares estimator (LSE)
with the algorithm of the R package EDR and the a penalized least squares estimator (PLSE),
using smoothing splines. The EDR package takes much longer than the other algorithms for
larger sample sizes, but is not superior to the other algorithms. In fact, an example is
given in BoxPlot_alpha_err_and_time_table1_n=5000.pdf, where the non-efficient simple
score estimator SSE gives a better performance, both in computing time and values.
I could not remove the text that is produced during the runs for the EDR package.

There exists an R package "simest" for the PLSE, but for simplicity I programmed the
procedure here from first principles, using the Reinsch algorithm as explained in the book
"Nonparametric Regression and Generalized Linear Models" of Green and Silverman (1994),
using systematically the band matrix structure to speed things up. The simple C++ code is
given in spline.cpp.

At the end of the runs of simulation_table1.R and simulation_table2.R the estimate of psi
computed in spline.cpp is shown. If one experiments with other sample sizes and
dimensions (by changing parameters at the beginning of the files
simulation_table1.R and simulation_table2.R), one might want to change the bandwidth or
the smoothing parameter mu in spline.cpp.

The computation of SSE, LSE and PLSE (the spline procedure) is started from alpha equal to
(1,0,...,0) (rather far from the actual value (alpha=1/sqrt{m})(1,...,1)), ESE is started
from the estimate given by SSE, and EDR does not need a starting value. The minimization
algorithm of Nelder-Mead is used instead of the Hooke-Jeeves method, used before.

For running the scripts, one needs to be able to run the package Rcpp. One also needs to
copy the other files of the directory, since most of them are used in the scripts.
Note that the estimates SSE, ESE and LSE expect data from a model where psi is monotone
increasing; otherwise erratic results can be expected. Also note that the efficient score
estimator can indeed give a significant improvement on SSE, when started at the result of
the SSE.
