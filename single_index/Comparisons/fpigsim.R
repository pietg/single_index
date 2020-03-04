#
##====================================================
## a fixed ponit iterative procedure for EFM for GSIM
##====================================================
#
## model     ---  E(y|x) = mu(g(beta*x))    
#               Var(y|x) = V(g(beta*x)) 
## arguments ---
#         xx --- predictors
#          y --- response 
#          r --- the deleted position, r in 1:p
#     Pconst --- the relaxation factor in Step 2, named by "M" 
#      index --- the initial parameter, default calue is (1,1,...,1)/sqrt(p)
#    mymodel --- designed for examples; "none" ( y follows normal distribution)
#                                       "binomial" ( y follows binomial distribution) 
#
##==============================================================
#
  fpigsim <- function (xx, y, r, Pconst, index, mymodel){
#----------------------------------------------------------------
  n <- dim(xx)[1]
  p <- dim(xx)[2]

  if (missing(Pconst))   {     
     Pconst <- 1/sqrt(p) }
  #if (missing(index))   {     
  #   index  <- rep(1,p)/sqrt(p) }
   if (mymodel=="binomial"){   
      # glmfit <- glm.fit(xx,y,family=binomial())
      # index0 <- glmfit$coefficients
      # index  <- index0/sqrt(index0%*%index0)  }
        index <- rep(1,p)/sqrt(p) }
   if (mymodel=="none") {
   index <- rep(1,p)/sqrt(p)}
   
#--------------------------------------------------------------- 
  para_old <- index   # the original value
  para_new <- rep(0,p)
   
   ep   <- 5e-5;                    # the  
 it_max <- 400;

    k   <- 1;

#---------------------------------------------------------------
   xxfit   <- array(0, dim=c(n,p))
   FF      <- rep(0, p)
   yfit    <- rep(0, n)
   deryfit <- rep(0, n)
   newx    <- rep(0, n)
#---------------------------------------------------------------
  while (k<=it_max){
   
   if (mymodel=="none"){    
         
        yfit0    <- smooth.lf(xx%*%para_old, y)
        yfit     <- yfit0$y
        deryfit0 <- smooth.lf(xx%*%para_old, y, deriv=1)
        deryfit  <- deryfit0$y;
       
#        ghat = npreg(tydat=y, txdat=xx%*%para_old, ckertype="epanechnikov", gradients=TRUE)
#        yfit = ghat$mean
#        deryfit = ghat$grad
       
       
        for (ii in 1:p){
            xxfit0       <- smooth.lf(xx%*%para_old, xx[ ,ii]) 
            xxfit[ , ii] <- xxfit0$y
            FF[ii]       <- sum((y-yfit)*deryfit*(xx[ ,ii] - xxfit[ , ii]))
        }    
     }

     if (mymodel=="binomial"|mymodel=="real-binomial"){   
        
        newx0 <- xx%*%para_old
        newx  <- newx0[,1]      
     
        yfit0    <- smooth.lf(newx, y, family="binomial", link="logit")
        yfit     <- yfit0$y
        deryfit0 <- smooth.lf(newx, y, family="binomial", link="logit", deriv=1)
        deryfit  <- deryfit0$y   
                 
       
        for (ii in 1:p){
            xxfit0       <- smooth.lf(newx, xx[ ,ii]) 
            xxfit[ , ii] <- xxfit0$y
        
            FF[ii]      <- sum((y-yfit)*deryfit*(xx[ ,ii] - xxfit[ , ii]));

           #  in0 <- which(y==0)
           #  in1 <- which(y==1)
                            

 
            # if (ii==r) {
            #         FF[ii]      <- sum(yfit[in0]*deryfit[in0]*(xx[in0, r] - xxfit[in0, r])) }
            # if (ii!=r) {
            #        FF[ii]      <-  sum(yfit[in0]*deryfit[in0]*(xx[in0, ii] - xxfit[in0, ii]))
            #                        -sum(yfit[in1]*deryfit[in1]*(xx[in1, ii] - xxfit[in1, ii]))
            #                        + para_old[ii]*sum((1-yfit[in1])*deryfit[in1]*(xx[in1, r] - xxfit[in1, r]))/para_old[r]  }

        }             
   
    }    
    
               

#----------------------------------------------------------------------------

        con      <- FF[r]/sqrt(crossprod(FF)); con <- con[1, 1]
        para_new <- para_old*Pconst/(Pconst+con) + FF*abs(con)/(Pconst+con)/sqrt(crossprod(FF))
        para_new <- para_new/sqrt(crossprod(para_new))
               
        norm     <- max(abs(para_new-para_old))
        if (norm<ep){ break }

        k <- k+1;     

        para_old <- para_new;
          
   } 
 
   list(root=para_new, it=k) 
}   
