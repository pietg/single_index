#
##====================================================
## a fixed point algorithm refined by fisher scoring method
##====================================================
#
## model     ---  E(y|x) = mu(g(beta*x))
#               Var(y|x) = V(g(beta*x))
## arguments ---
#         xx --- predictors
#          y --- response
#          r --- the deleted position, r in 1:p
#    mymodel --- designed for examples; "none" ( y follows normal distribution)
#                                       "binomial" ( y follows binomial distribution)
#
##==============================================================
#
  fisher <- function (xx, y, r, mymodel){

#----------------------------------------------------------------
  n <- dim(xx)[1]
  p <- dim(xx)[2]

#---------------------------------------------------------------

  para_old <- rep(0, p)

  if (mymodel=="none"){

       para_old0  <- fpigsim(xx, y, r, Pconst=2.5*(p<=10)+0.5*(p>10), index=rep(1,dim(xx)[2])/sqrt(dim(xx)[2]), mymodel="none")
       para_old   <- para_old0$root
        }

       # Here "Pconst" plays the role of "M" in paper.
       # We fix "Pconst" for each setting in our examples. Though it maybe lose some accuracy, it saves much computation time.


  if (mymodel=="binomial"){
      #-------------------- Design C -----------------
      #para_old0  <- fpigsim(xx, y, r, Pconst=0.46*(p>10) + 2.5*(p<=10), index=rep(1,p)/sqrt(p), mymodel="binomial")
      #para_old   <- para_old0$root
      #-------------------- Design D -----------------
      para_old0  <- fpigsim(xx, y, r, Pconst=1.1*(p>10) + 5.5*(p<=10), index=rep(1,p)/sqrt(p), mymodel="binomial")
      para_old   <- para_old0$root

       }

      # Here "Pconst" plays the role of "M" in paper.
      # We fix "Pconst" for each setting in our examples. Though it maybe lose some accuracy, it saves much computation time.


  if (mymodel=="real-binomial"){

     para_old0  <- fpigsim(xx, y, r, Pconst=1.3/sqrt(p), index=rep(1,p)/sqrt(p), mymodel="binomial")

     para_old   <- para_old0$root }

#--------------------------------------------

  ep        <- 1e-6
  it_max    <- 800*(mymodel=="none") + 0*(mymodel!="none")
  #it_max <- 800
  k         <- 1

#---------------------------------------------------------------

   xxfit     <- array(0, dim=c(n,p));
   yfit      <- rep(0, n);
   deryfit   <- rep(0, n);
   para_new  <- rep(0, p)
   Jac       <- array(0, dim=c(p, p-1))
#---------------------------------------------------------------
  while (k<=it_max){

   if (mymodel=="none"){

       #--------------------------------------------------
        yfit0    <- smooth.lf(xx%*%para_old, y)
        yfit     <- yfit0$y
        deryfit0 <- smooth.lf(xx%*%para_old, y, deriv=1)
        deryfit  <- deryfit0$y;
        for (ii in 1:p){
            xxfit0       <- smooth.lf(xx%*%para_old, xx[ ,ii])
            xxfit[ , ii] <- xxfit0$y  }
       #---------------------------------------------------

        Jac <- rbind(-para_old[-r]/sqrt(1-crossprod(para_old[-r])), diag(p-1))

        FisherM <- t(Jac)%*%(t(deryfit*(xx - xxfit))%*%(deryfit*(xx - xxfit))/n)%*%Jac

        EQ     <-  t(Jac)%*%colMeans(as.vector((y-yfit))*deryfit*(xx - xxfit))
        
     }

     if (mymodel=="binomial"|mymodel=="real-binomial"){

        #--------------------------------------------------
        yfit0    <- smooth.lf(xx%*%para_old, y, family="binomial", link="logit")
        yfit     <- yfit0$y
        deryfit0 <- smooth.lf(xx%*%para_old, y, family="binomial", link="logit", deriv=1)
        deryfit  <- deryfit0$y;
        for (ii in 1:p){
            xxfit0       <- smooth.lf(xx%*%para_old, xx[ ,ii])
            xxfit[ , ii] <- xxfit0$y  }
       #---------------------------------------------------

        Jac <- rbind(-para_old[-r]/sqrt(1-crossprod(para_old[-r])), diag(p-1))

        FisherM <- t(Jac)%*%(t(deryfit*(xx - xxfit))%*%(deryfit*(xx - xxfit))/n)%*%Jac

        EQ     <-  t(Jac)%*%colMeans((y-yfit)*deryfit*(xx - xxfit))

      }



#----------------------------------------------------------------------------

     if (kappa(FisherM)<5e2){

        para_new[-r] <- para_old[-r] + ginv(FisherM)%*%EQ
        para_new[r]  <- para_old[r]
        para_new <- para_new/sqrt(crossprod(para_new)) }

     norm     <- max(abs(para_new-para_old))

     if (kappa(FisherM)>5e2|norm<ep){ break }

     k <- k+1;

     para_old <- para_new;

   }

  # para_new[1:2] <- para_old[1:2]/sqrt(para_old[1:2]%*%para_old[1:2])
  # para_new[3:p] <- (para_old[3:p])^2
  # para_new <- para_new/sqrt(para_new%*%para_new)

    para_new <- para_old

   list(root=para_new, root_fp=para_old0$root, it=k)
}
