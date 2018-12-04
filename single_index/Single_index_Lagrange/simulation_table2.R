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
sourceCpp("MRE.cpp")
sourceCpp("linear.cpp")
sourceCpp("spline.cpp")

# NumIt is the number of replications
  NumIt = 100
  
# n is the number of observations
  n = 200
  
# m is the dimension of the covariate
  m= 15

 # sigma is the variance of the error
  sigma = 1
 
# a0 is alpha_0, the finite dimensional parameter
  a0 = c(rep(1,m))/sqrt(m)
  
  Sigma2 <- m*diag(m)+matrix(rep(-1,times=m^2),nrow=m,ncol=m)
  Sigma2 <- (1/(m*sqrt(27)))*Sigma2
  #Sigma2 <- (1/(m*sqrt(34)))*Sigma2
  
  timeMat <- NULL
  normMat <- matrix(0, nrow= NumIt, ncol= 8)
  colnames(normMat) <- c("EDR","SSE","ESE","LSE","MRE","linear","PLSE","PLSE2")

for (j in 1: NumIt){
  
  sim = 1001+j
   
  set.seed(sim)
  
  print(j)	
	X = matrix(rnorm(m*n,0,sigma),n,m, byrow = FALSE)
	z=X%*%a0
	y= z^3 + rnorm(n,0,sigma)

	# EDR estimate proposed by Hristache et al.
	starter_edr = proc.time()
	EDR <- edr(X,y,method = "HJPS")
	edr_hat = -summary(EDR)$Rhat[1,]
	time_edr = (proc.time() -starter_edr)[3]
	
	# LSE
	starter_lse = proc.time()
	LSE <- ComputeLSE(X,y,20,m)
	lse_hat = LSE$alpha
	time_lse = (proc.time() -starter_lse)[3]
	
	# linear
	starter_linear = proc.time()
	linear <- Compute_linear(X,y,lse_hat,m)
	linear_hat = linear$alpha
	time_linear = (proc.time() -starter_linear)[3]
	
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
	
	# MRE
	starter_mre = proc.time()
	MRE <- ComputeMRE(X,y,lse_hat,m)
	mre_hat = MRE$alpha
	time_mre = (proc.time() -starter_mre)[3]
	
	# PLSE
	starter_PLSE = proc.time()
	PLSE <- Compute_spline(X,y,lse_hat,m)
	PLSE_hat <- PLSE$alpha
  	time_PLSE = (proc.time() -starter_PLSE)[3]
  
# PLSE from R package simest, with LSE as starting point
    starter_PLSE2 = proc.time()
    PLSE2 <- sim.est(X,y, method = "smooth.pen", beta.init = lse_hat, lambda = .1, maxit = 1000)
    PLSE_hat2 <- PLSE2$beta
    time_PLSE2 = (proc.time() -starter_PLSE2)[3]
	
	  # asymptotic distribution
	starter_asymp= proc.time()
  	x1 <- rnorm(m,0,1)
  	x2<-Sigma2%*%x1
  	asymp<-norm(x2,"2")/sqrt(m)
  	time_asymp = (proc.time() -starter_asymp)[3]
  	
  write(edr_hat,file = "EDR.txt", ncol =m, append = TRUE)
  write(ese_hat,file = "ESE.txt", ncol =m, append = TRUE)
  write(sse_hat,file = "SSE.txt", ncol =m, append = TRUE)
  write(lse_hat,file = "LSE.txt", ncol =m, append = TRUE)
  write(linear_hat,file = "linear.txt", ncol =m, append = TRUE)
  write(mre_hat,file = "MRE.txt", ncol =m, append = TRUE)
  write(PLSE_hat,file = "PLSE.txt", ncol =m, append = TRUE)
  write(PLSE_hat2,file = "PLSE2.txt", ncol =m, append = TRUE)
   
  	normMat[j,]  = c(sqrt(n)*norm((edr_hat- a0), "2")/sqrt(m),
sqrt(n)*norm((sse_hat- a0), "2")/sqrt(m),
sqrt(n)*norm((ese_hat- a0), "2")/sqrt(m),
sqrt(n)*norm((lse_hat- a0), "2")/sqrt(m),
sqrt(n)*norm((mre_hat-a0), "2")/sqrt(m),
sqrt(n)*norm((linear_hat-a0), "2")/sqrt(m),
sqrt(n)*norm((PLSE_hat-a0), "2")/sqrt(m),
sqrt(n)*norm((PLSE_hat2-a0), "2")/sqrt(m))
	timeMat<-rbind(timeMat,c(time_edr,time_sse,time_ese,time_lse,time_mre,time_linear,time_PLSE,time_PLSE2))
}

colnames(timeMat) <- c("EDR","SSE","ESE","LSE","MRE","linear","PLSE","PLSE2")
pdf("BoxPlot_alpha_err_and_time.pdf")
boxplot(normMat, main= "Boxplot of sqrt(n/d)||alpha_hat-alpha_0||_2", las=2)
#boxplot(normMat)
#boxplot(timeMat, main="Run Times", las=2) 
dev.off()

	A <- PLSE$psi
	B <- PLSE$derivative
	C <- PLSE$data

    x1<-A[,1]
    y1<-A[,2]
    x2<-B[,1]
    y2<-B[,2]
    x<-C[,1]
   	y<-C[,2]

    f <- function(x) {x^3}
   	x0<-seq(min(x1),max(x1),by=0.01)
    y0<-f(x0)
    plot(c(-10000,-10000),xlim=c(min(x1),max(x1)), ylim=c(min(y1,y0),max(y1,y0)), main= "",ylab="",xlab="",bty="n",las=1)
    lines(x1,y1,col="red",lwd=2)
    lines(x0,y0,lwd=2,col="blue",lty=2)
    points(x,y,pch = 21)
    
    g <- function(x) {3*x^2}
    x0<-seq(min(x2),max(x2),by=0.1)
    y0<-g(x0)
    plot(c(-100,-100),xlim=c(min(x2),max(x2)), ylim=c(min(y2,y0),max(y2,y0)), main= "",ylab="",xlab="",bty="n",las=1)
    lines(x2,y2,lwd=2)
    lines(x0,y0,lwd=2,col="red",lty=2)
    
    
    
    B<-read.table("EDR.txt")
    mean(B[,1])
    mean(B[,2])
    mean(B[,3])
	n*var(B)
	B<-read.table("SSE.txt")
	mean(B[,1])
    mean(B[,2])
    mean(B[,3])
	n*var(B)
	B<-read.table("ESE.txt")
	mean(B[,1])
    mean(B[,2])
    mean(B[,3])
	n*var(B)
	B<-read.table("LSE.txt")
	mean(B[,1])
    mean(B[,2])
    mean(B[,3])
	n*var(B)
	B<-read.table("linear.txt")
	mean(B[,1])
    mean(B[,2])
    mean(B[,3])
	n*var(B)
	B<-read.table("MRE.txt")
	mean(B[,1])
    mean(B[,2])
    mean(B[,3])
	n*var(B)
	B<-read.table("PLSE.txt")
	mean(B[,1])
    mean(B[,2])
    mean(B[,3])
	n*var(B)
	B<-read.table("PLSE2.txt")
	mean(B[,1])
    mean(B[,2])
    mean(B[,3])
	n*var(B)
    
    