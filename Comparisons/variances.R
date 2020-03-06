########################################################
####   MONOTONE REGRESSION        ######
########################################################
n=500

B<-read.table("EDR.txt")
colMeans(B)
n*cov(B)

B<-read.table("LSE.txt")
colMeans(B)
n*cov(B)

B<-read.table("ESE.txt")
colMeans(B)
n*cov(B)

B<-read.table("SSE.txt")
colMeans(B)
n*cov(B)

B<-read.table("spline.txt")
colMeans(B)
n*cov(B)

B<-read.table("MAVE.txt")
colMeans(B)
n*cov(B)

B<-read.table("EFM.txt")
colMeans(B)
n*cov(B)


