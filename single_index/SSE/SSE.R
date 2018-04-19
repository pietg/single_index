	library(Rcpp)
	sourceCpp("SSE.cpp")
  output <- Compute_SSE()
   
	A <- output$SSE
	B <- output$means
	C <- output$covariance_matrix
	
	B
	C
