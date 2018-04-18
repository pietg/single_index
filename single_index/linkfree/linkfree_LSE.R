	library(Rcpp)
	sourceCpp("linkfree_LSE.cpp")
  output <- ComputeLinkFree_LSE()
   
	A <- output$LFLSE
	B <- output$means
	C <- output$covariance_matrix
	
	B
	C
