########################################################
####   MONOTONE REGRESSION        ######
########################################################
rm(list=ls())
source('fpigsim.R')
source('fisher.R')
library(locfit)
library(Rcpp)
library(MAVE)
library(MASS)
sourceCpp("SSE.cpp")
sourceCpp("ESE.cpp")
sourceCpp("LSE.cpp")
sourceCpp("spline.cpp")

  NumIt = 100
  n = 100
  m= 2
  sigma = 1
  a0 = c(1/sqrt(2),1/sqrt(2))
  
  timeMat <- NULL
   normMat <- matrix(0, nrow= NumIt, ncol= 6)
  colnames(normMat) <- c("SSE","ESE","LSE","spline","EFM","MAVE")


for (j in 1: NumIt){
  
  sim = 1+j
   
  set.seed(sim)
  
  print(j)
	
	X = matrix(runif(m*n,0,1),n,m, byrow = FALSE)
	z=X%*%a0
	y=(z)^3 + rnorm(n,0,sigma)
	#y=rbinom(n,10,exp(z)/(1+exp(z)))

	# EDR estimate proposed by Hristache et al.
	#starter_edr = proc.time()
	#EDR <- edr(X,y,method = "HJPS")
	#edr_hat = -summary(EDR)$Rhat[1,] 
	#time_edr = (proc.time() -starter_edr)[3]
	
	# LSE
	starter_lse = proc.time()
	LSE <- ComputeLSE(X,y)
	lse_hat = LSE$alpha
	time_lse = (proc.time() -starter_lse)[3]
	
	# SSE
	starter_sse = proc.time()
	SSE <- ComputeSSE(X,y)
	sse_hat = SSE$alpha
	time_sse = (proc.time() -starter_sse)[3]
	
	# ESE
	starter_ese = proc.time()
	ESE <- ComputeESE(X,y)
	ese_hat = ESE$alpha
	time_ese = (proc.time() -starter_ese)[3]
		
	# spline
	starter_spline = proc.time()	
	spline <-Compute_spline(X,y)
	spline_hat = spline$alpha
	time_spline = (proc.time() -starter_spline)[3]

	# MAVE
	starter_MAVE = proc.time()
	MAVE <- mave.compute(X,y,method = 'meanmave')
	MAVE_hat = MAVE$dir[[1]]
	if (MAVE_hat[1]<0)
	{
		MAVE_hat[1]=-MAVE_hat[1]
		MAVE_hat[2]=-MAVE_hat[2]
	}
	time_MAVE = (proc.time() -starter_MAVE)[3]
	
# EFM
starter_EFM = proc.time()
EFM <- fisher(X,y,1,mymodel="none")
EFM_hat <- EFM$root
time_EFM = (proc.time()-starter_EFM)[3]

  write(ese_hat,file = "ESE.txt", ncol =m, append = TRUE)
  write(sse_hat,file = "SSE.txt", ncol =m, append = TRUE)
  write(lse_hat,file = "LSE.txt", ncol =m, append = TRUE)
  write(MAVE_hat,file = "MAVE.txt", ncol =m, append = TRUE)
  write(EFM_hat,file = "EFM.txt", ncol =m, append = TRUE)
  write(spline_hat,file = "spline.txt", ncol =m, append = TRUE)
  
  normMat[j,]  = c(norm(sse_hat- a0, "2"),norm(ese_hat- a0, "2"),norm(lse_hat-a0, "2"),norm(spline_hat-a0, "2"),norm(EFM_hat-a0, "2"),norm(MAVE_hat-a0, "2"))
timeMat<-rbind(timeMat,c(time_sse,time_ese,time_lse,time_spline,time_MAVE))
}

colnames(timeMat) <- c("SSE","ESE","LSE","spline","MAVE")
pdf("BoxPlot_alpha_err.pdf")
boxplot(normMat,las=1)
boxplot(timeMat, main="Run Times", las=1) 
dev.off()

	A <- ESE$psi
	C <- ESE$data

   x1<-A[,1]
    y1<-A[,2]
    x<-C[,1]
   	y<-C[,2]

    f <- function(x) {x^3}
    #f <- function(x) {10*exp(x)/(1+exp(x))}
    x0<-seq(min(x1,x),max(x1,x),by=0.01)
    y0<-f(x0)
    plot(c(-10000,-10000),xlim=c(min(x1,x),max(x1,x)), ylim=c(min(y,y0,y1),max(y,y0,y1)), main= "",ylab="",xlab="",bty="n",las=1)
    lines(x1,y1,col="blue",lwd=2,type='s')
    lines(x0,y0,lwd=2,col="red",lty=2)
    points(x,y,pch = 21)
