########################################################
####   MONOTONE REGRESSION        ######
########################################################
rm(list=ls())
library(Rcpp)
library(EDR)
source("NewSimEst.R")
sourceCpp("SSE.cpp")
sourceCpp("ESE.cpp")
sourceCpp("LSE.cpp")

  NumIt = 100
  n = 100
  m= 25
  sigma = 1
  a0 = c(rep(1,m))/sqrt(m)
  
  Sigma2 <- m*diag(m)+matrix(rep(-1,times=m^2),nrow=m,ncol=m)
  Sigma2 <- (1/(m*sqrt(27)))*Sigma2
  
  timeMat <- NULL
  normMat <- matrix(0, nrow= NumIt, ncol= 6)
  colnames(normMat) <- c("EDR", "SSE","ESE","LSE","PlSE","asymp")

for (j in 1: NumIt){
  
  sim = 101+j
   
  set.seed(sim)
  
  print(j)
	
	X = matrix(rnorm(m*n,0,sigma),n,m, byrow = FALSE)
	z=X%*%a0
	y=(z)^3 + rnorm(n,0,sigma)

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
	SSE <- ComputeSSE(X,y,a0,m)
	sse_hat = SSE$alpha
	time_sse = (proc.time() -starter_sse)[3]
	
	# ESE
	starter_ese = proc.time()
	ESE <- ComputeESE(X,y,a0,m)
	ese_hat = ESE$alpha
	time_ese = (proc.time() -starter_ese)[3]

  # PLSE with a0 as starting point
  starter_PLSE= proc.time()
  PLSE_hat <- Alter_Min_Simest(X, y, beta_init =a0, lambda = 0.1, method ="smooth.pen", nmax= 100, maxit=100)
  time_PLSE = (proc.time() -starter_PLSE)[3]
	
	  # asymptotic distribution
	starter_asymp= proc.time()
  	x1 <- rnorm(m,0,1)
  	x2<-Sigma2%*%x1
  	asymp<-norm(x2,"2")/sqrt(m)
  	time_asymp = (proc.time() -starter_asymp)[3]

  	normMat[j,]  = c(sqrt(n)*norm((edr_hat- a0), "2")/sqrt(m),sqrt(n)*norm((sse_hat- a0), "2")/sqrt(m),sqrt(n)*norm((ese_hat- a0), "2")/sqrt(m),sqrt(n)*norm((lse_hat-a0), "2")/sqrt(m),sqrt(n)*norm((PLSE_hat-a0), "2")/sqrt(m),asymp)
	timeMat<-rbind(timeMat,c(time_edr,time_sse,time_ese,time_lse,time_PLSE,time_asymp))
}

colnames(timeMat) <- c("EDR","SSE","ESE","LSE","PLSE","asymp")
pdf("BoxPlot_alpha_err_and_time.pdf")
boxplot(normMat, main= "Boxplot of sqrt(n/d)||alpha_hat-alpha_0||_2", las=2)
boxplot(timeMat, main="Run Times", las=2) 
dev.off()

	A <- ESE$psi
	B <- ESE$data

    x1<-A[,1]
    y1<-A[,2]
    x<-B[,1]
   	y<-B[,2]

    f <- function(x) {x^3}
    x0<-seq(min(x1),max(x1),by=0.01)
    y0<-f(x0)
    plot(c(-10000,-10000),xlim=c(min(x1),max(x1)), ylim=c(min(y,y0,y1),max(y,y0,y1)), main= "",ylab="",xlab="",bty="n",las=1)
    lines(x1,y1,col="blue",lwd=2,type="s")
    lines(x0,y0,lwd=2,col="red",lty=2)
    points(x,y,pch = 21)