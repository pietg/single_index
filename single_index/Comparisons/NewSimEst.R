## foo computes the objective value for any given value
## of the index parameter beta.
## method can be one of the following values
## 	(1) cvx.lse
## 	(2) cvx.lip
## 	(3) cvx.pen
## 	(4) smooth.pen

library("simest")
## function to compute the least squares loss for a given index parameter.
## and a few additional calculations.
foo <- function(beta, xx, yy, method, lambda, L, flag = 0){
	xx_beta <- xx%*%beta
	A <- cbind(xx_beta, yy, xx)
	A <- A[order(A[,1]),]
	xx_beta_new <- A[,1]
	yy_new <- A[,2]
	xx_new <- A[,-c(1,2)]
	command <- paste0(method, ".reg(xx_beta_new, yy_new, lambda = lambda, L = L)")
	tmp <- eval(parse(text = command))
	if(flag == 0){
		G <- NULL
		R1 <- NULL
		R2 <- NULL
		beta_tau_pos <- NULL
		beta_tau_neg <- NULL
	}
	if(flag == 1){
		G <- colMeans(tmp$residuals*tmp$deriv*xx_new)
		Beta_G <- sum(beta*G)
		GbyBeta <- G[1]/beta[1]
		bb <- Beta_G - GbyBeta
		ac <- norm(G, "2")^2 - Beta_G^2
		R1 <- -(bb - sqrt(bb^2 + ac))/(0.5*ac)
		R2 <- -(bb + sqrt(bb^2 + ac))/(0.5*ac)
		beta_tau_pos <- function(alpha){
			tau <- alpha*R2
			(beta*(1 - 0.25*ac*tau^2 + tau*Beta_G) - tau*G)/(1 + 0.25*ac*tau^2)
		}
		beta_tau_neg <- function(alpha){
			tau <- alpha*R1
			(beta*(1 - 0.25*ac*tau^2 + tau*Beta_G) - tau*G)/(1 + 0.25*ac*tau^2)
		}
	}
	return(list(fval = sum(tmp$residuals^2), gradient = G, R1 = R1, R2 = R2, beta_tau_pos = beta_tau_pos, beta_tau_neg = beta_tau_neg))
}

## finding tau that reduces the function value; no minimization for tau.
find_tau <- function(beta, xx, yy, method, lambda, L, nmax = 10, tol = 1e-05){
	if(beta[1] < tol) return(beta)
	else{
		run_beta <- foo(beta, xx, yy, method, lambda, L, flag = 1)
		R1 <- run_beta$R1; R2 <- run_beta$R2
		if(abs(R2 - R1) < tol) return(beta)
		fval0 <- run_beta$fval
		BETA_TAU_POS <- function(alpha) run_beta$beta_tau_pos(alpha)
		FOO_pos <- function(alpha){
			foo(BETA_TAU_POS(alpha), xx, yy, method, lambda, L)$fval
		}
		alpha_start <- 0.99
		for(i in 1:nmax){
			if(FOO_pos(alpha_start) < fval0){
				return(BETA_TAU_POS(alpha_start))
			} else {
				alpha_start <- alpha_start/2
			}
		}
		BETA_TAU_NEG <- function(alpha) run_beta$beta_tau_neg(alpha)
		FOO_neg <- function(alpha){
			foo(BETA_TAU_NEG(alpha), xx, yy, method, lambda, L)$fval
		}
		alpha_start <- 0.99
		for(i in 1:nmax){
			if(FOO_neg(alpha_start) < fval0){
				return(BETA_TAU_NEG(alpha_start))
			} else {
				alpha_start <- alpha_start/2
			}
		}
	}
	return(beta)
}

Alter_Min_Simest <- function(xx, yy, beta_init, method, lambda, L, nmax = 10, maxit = 100, beta_tol = 1e-05, beta_pos_tol = 1e-05){
	beta_start <- beta_end <- beta_init
	flag <- 0; i <- 1
	while(flag == 0 && i <= maxit){
		beta_end <- find_tau(beta_start, xx, yy, method = method, lambda, L, nmax = nmax, tol = beta_pos_tol)
		if(norm(beta_start - beta_end, "2") < beta_tol){
			flag <- 1
		} else {
			beta_start <- beta_end
		}
		i <- i + 1 
	}
	if(flag == 0){
		print(paste0("Algorithm terminated by reaching", maxit, "alternating steps."))
	}
	#warning("Algorithm terminated since it reached maxit")
	return(beta_end)
}

set.seed(1)
n <- 100; d <- 3
beta0 <- rep(1, d)/sqrt(d)
m0 <- function(t){ t^2 }
eps <- rnorm(n, 0, 0.2)
xx <- matrix(runif(n*d, -1, 1), ncol = d)
yy <- m0(xx%*%beta0) + eps
beta_init <- c(1,0,0)
Alter_Min_Simest(xx, yy, beta_init, method = "cvx.pen", lambda = 0.01, nmax = 101, maxit=13)
## somehow nmax = 100 always returns 0 to me.
## why? what is the mistake in this code?