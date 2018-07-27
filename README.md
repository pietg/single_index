# single_index

This repository gives Rcpp scripts for simulations with the simple score estimator and a
link-free estimate, using the methods of section 4.2 and 5 of "Score estimation in the
monotone single index model" by Fadoua Balabdaoui, Piet Groeneboom and Kim Hendrickx,
https://arxiv.org/abs/1712.05593

Then R scripts simulation_table1.R (for Table 1 in the paper) and simulation_table2.R
(for Table 2 in the paper) in the directory ``Comparisons'' compare the simple score
estimator (SSE), the efficient score estimator (ESE) and the least squares estimator (LSE)
with the algorithm of the R package EDR and the penalized least squares estimator (PLSE)
of the R package simest. The EDR package takes much longer than the other algorithms for
larger sample sizes, but is not superior to the other algorithms. In fact, an example is
given in BoxPlot_alpha_err_and_time_table1_n=5000.pdf, where the non-efficient simple
score estimator SSE gives a better performance, both in computing time and values.
I could not remove the text that is produced during the runs for the EDR package.

At the end of the runs of simulation_table1.R and simulation_table2.R the estimate of psi
and its derivative psi' computed in ESE is shown. If one experiments with other sample
sizes and dimensions (by changing parameters at the beginning of the files
simulation_table1.R and simulation_table2.R) the bandwidth for estimating psi' has to be
adapted in ESE.cpp.

The computation of SSE and LSE is started from alpha=(1,0,...,0) (rather far from the
actual value (alpha=1/sqrt{m})(1,...,1)), the other algorithms except EDR are started from
the estimate given by SSE, and EDR does not need a starting value. The minimization
algorithm of Nelder-Mead is used instead of the Hooke-Jeeves method, used before.
For running the scripts, one needs to be able to run the package Rcpp. One also needs to
copy the other files of the directory, since most of them are used in the scripts.
Note that the estimates SSE,ESE and LSE expect data from a model where psi is monotone
increasing; otherwise erratic results can be expected. Also note that the efficient score
estimator can indeed give a significant improvement on SSE, when started at the result of
the SSE.
