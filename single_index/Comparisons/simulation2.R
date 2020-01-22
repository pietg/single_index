########################################################
####   MONOTONE REGRESSION        ######
########################################################
rm(list=ls())
library(Rcpp)
library(EDR)
sourceCpp("SSE.cpp")
sourceCpp("ESE.cpp")
sourceCpp("LSE.cpp")
sourceCpp("spline.cpp")

  NumIt = 1000
  n = 100
  m= 2
  sigma = 1
  
  normMat <- matrix(0, nrow= NumIt, ncol= 5)
  colnames(normMat) <- c("EDR", "SSE","ESE","LSE","spline")


for (j in 1: NumIt){
  
  sim = 2001+j
   
  set.seed(sim)
  
  print(j)

  a0 = c(1/sqrt(2),1/sqrt(2))
	
	X = matrix(runif(m*n,0,1),n,m, byrow = FALSE)
	#X = matrix(rnorm(m*n,0,sigma),n,m, byrow = FALSE)
	z=X%*%a0
	#y=(z)^3 + rnorm(n,0,sigma)
	y=rbinom(n,10,exp(z)/(1+exp(z)))

	# EDR estimate proposed by Hristache et al.
	starter_edr = proc.time()
	EDR <- edr(X,y,method = "HJPS")
	edr_hat = -summary(EDR)$Rhat[1,]
	time_edr = (proc.time() -starter_edr)[3]
	
	# LSE
	starter_lse = proc.time()
	LSE <- ComputeLSE(X,y,a0,m)
	lse_hat = LSE$alpha
	time_lse = (proc.time() -starter_lse)[3]
	
	# SSE
	starter_sse = proc.time()
	SSE <- ComputeSSE(X,y,lse_hat,m)
	sse_hat = SSE$alpha
	time_sse = (proc.time() -starter_sse)[3]
	
	# ESE
	starter_ese = proc.time()
	ESE <- ComputeESE(X,y,lse_hat,m)
	ese_hat = ESE$alpha
	time_ese = (proc.time() -starter_ese)[3]
		
	# spline
	starter_spline = proc.time()	
	spline <-Compute_spline(X,y,lse_hat,m)
	spline_hat = spline$alpha
	time_spline = (proc.time() -starter_spline)[3]

  write(edr_hat,file = "EDR.txt", ncol =m, append = TRUE)
  write(ese_hat,file = "ESE.txt", ncol =m, append = TRUE)
  write(sse_hat,file = "SSE.txt", ncol =m, append = TRUE)
  write(lse_hat,file = "LSE.txt", ncol =m, append = TRUE)
  write(spline_hat,file = "spline.txt", ncol =m, append = TRUE)
  
  normMat[j,]  = c(norm(edr_hat- a0, "2"),norm(sse_hat- a0, "2"),norm(ese_hat- a0, "2"),norm(lse_hat-a0, "2"),norm(spline_hat-a0, "2"))
}

#colnames(timeMat) <- c("EDR","SSE","ESE","LSE","spline")
pdf("BoxPlot_alpha_err.pdf")
boxplot(normMat,las=1)
dev.off()

