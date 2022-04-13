###################################################################################################################
######   Flexible joint modeling of mean and dispersion for the directional tuning of neuronal spike counts  ######
######                     								                                                                   ######
######                    María Alonso-Pena,	    Irène Gijbels    and   Rosa M. Crujeiras                   ######
######                                                                          								             ######
######                     Universidade de         KU Leuven		          Universidade de                    ######
######   		                  Santiago de 				                          Santiago de                      ######                                  
######  		                  Compostela                                    Compostela                       ######
######                     								                                                                   ######
###################################################################################################################

## Supplementary R code for the paper 
## "Flexible joint modeling of mean and dispersion for the directional tuning of neuronal spike counts" 

#################
library(circular)
library(Rcpp)
#################

# NONPARAMETRIC DOUBLE POISSON MEAN-AND-DISPERSION ESTIMATION FOR A CIRCULAR COVARIATE

# The function kern.dp.circ performs the nonparametric joint estimation of the mean
# and dispersion functions when the covariate is a circular variable and the response
# follows a conditional double Poisson distribution. The estimation is carried out
# with a local linear type sine-polynomial estimator (p=1 in equations (13) and (14) 
# of Alonso-Pena et al. (2022)). 

# The function uses the internal code in "DoublePoisson.cpp", for which it is necessary
# to load the Rcpp package. The arguments are the following:

## t: vector with points where the regression functions are estimated. Units must be radians.
## x: vector of data for the independent variable. Units must be radians.
## y: vector of data for the dependent variable. The must represent counts.
## kappa: smoothing parameter to be used in the estimation of the mean.
## nu: smoothing parameter to be used in the estimation of the dispersion.

# It returns a list with two components: 
# A vector with the estimations of g() at t
# A vector with the estimations of psi() at t

sourceCpp("DoublePoisson.cpp")

kern.dp.circ<-function(t,x,y,kappa,nu){
  
  if (!is.numeric(x))
    stop("argument 'x' must be numeric")
  if (!is.numeric(y))
    stop("argument 'y' must be numeric")
  if (!is.numeric(t))
    stop("argument 't' must be numeric")
  if (length(x) != length(y))
    stop("'x' and 'y' must have the same number of observations")
  if (!is.numeric(kappa))
    stop("argument 'kappa' must be numeric")
  if (kappa<=0)
    stop("argument 'kappa' must be greater than 0")
  if (!is.numeric(nu))
    stop("argument 'nu' must be numeric")
  if (nu<=0)
    stop("argument 'nu' must be greater than 0")
  
  startv<-c(log(mean(y)),0,0,0)
  startv2<-c(0,0)
  
  est<-my_fun_DoublePois(t,x,y,startv,startv2,kappa,nu)
  
  return(est)
  
}




# AUTOMATIC SELECTION OF THE SMOOTHING PARAMETERS VIA THE TWO-STEP CROSS-VALIDATION


# AUXILIARY FUNCTIONS: pois_p1 and loglik_cv

# Efficient estimation p=1 JUST MEAN

pois_p1<-function(xseq,x,y,kappa,startv){
  
  k<-length(xseq)
  beta0<-rep(startv[1],k)
  beta1<-rep(startv[2],k)
  
  x_t <- outer(x, xseq, "-")
  kxt <- exp(kappa * cos(x_t))/(besselI(kappa,0)*2*pi)
  s <- sin(x_t)
  ts <- t(beta0+beta1*t(s))
  ets <- exp(ts)
  bn0 <- apply(kxt * ets, 2, sum)
  bn1 <- apply(kxt * ets * s, 2, sum)
  bn2 <- apply(kxt * ets * s^2, 2, sum)
  bjti0 <-  t(kxt)*(bn2 - t(s) * bn1)
  #L0 <- t(t(bjti0)/(bn2*bn0-bn1^2))
  
  
  bjti1 <-  t(kxt)*(t(s)*bn0 - bn1)
  #L1 <- t(t(bjti1)/(bn2*bn0-bn1^2))
  
  K <- (y-ets)
  
  #beta0<-beta0 + diag(bjti0%*%K/(bn2*bn0-bn1^2))
  #beta1<-beta1 + diag(bjti1%*%K/(bn2*bn0-bn1^2))
  
  beta0<-beta0 + colSums(t(bjti0) * K)/(bn2*bn0-bn1^2)
  beta1<-beta1 + colSums(t(bjti1) * K)/(bn2*bn0-bn1^2)
  
  res<-5000
  tol<-0.001
  b<-0
  while(res>tol){
    
    
    beta0_old<-beta0
    beta1_old<-beta1
    
    ts <- t(beta0+beta1*t(s))
    ets <- exp(ts)
    K <- (y-ets)
    bn0 <- apply(kxt * ets, 2, sum)
    bn1 <- apply(kxt * ets * s, 2, sum)
    bn2 <- apply(kxt * ets * s^2, 2, sum)
    bjti0 <- t(kxt)*(bn2 - t(s) * bn1)
    bjti1 <- t(kxt)*(t(s)*bn0 - bn1)
    
    
    #beta0<-beta0 + diag(bjti0%*%K/(bn2*bn0-bn1^2))
    #beta1<-beta1 + diag(bjti1%*%K/(bn2*bn0-bn1^2))
    
    beta0<-beta0 + colSums(t(bjti0) * K)/(bn2*bn0-bn1^2)
    beta1<-beta1 + colSums(t(bjti1) * K)/(bn2*bn0-bn1^2)
    
    res<-max(c(abs(beta0-beta0_old),abs(beta1-beta1_old)))
    b<-b+1
    #print(b)
  }
  return(list(beta0,beta1))
}



#  cross-validation for p=1 (through log-likelihood) JUST MEAN

loglik_cv <- function(x,y,kappa,startv) {	
  n<-length(x)
  error <- numeric(n)
  for (j in 1:n) {
    
    fitj<-try(pois_p1(x[j],x[-j],y[-j],kappa,startv)[[1]],silent=TRUE)
    error[j] <- try(y[j]*fitj - exp(fitj),silent=TRUE)
  }
  me<-try(mean(error),silent=TRUE)
  if(class(me)!="try-error"){
    return(me)
  }
  
}


# TWO-STEP CROSS-VALIDATION 

# The function cv_kappa_nu performs the two-step cross-validation method described
# in Section 5.1 of Alonso-Pena et al. (2022). The arguments are the following:

## x: vector of data for the independent variable. Units must be radians.
## y: vector of data for the dependent variable. The must represent counts.
## lowerk: minimum value for the search of the smoothing parameter in mean estimation.
## upperk: maximum value for the search of the smoothing parameter in mean estimation.
## lowerknu: minimum value for the search of the smoothing parameter in dispersion estimation.
## upperknu: maximum value for the search of the smoothing parameter in dispersion estimation.

cv_kappa_nu<-function(x,y,lowerk=0.03,upperk=50,lowerknu=0.03,upperknu=7){
  
  if (!is.numeric(x))
    stop("argument 'x' must be numeric")
  if (!is.numeric(y))
    stop("argument 'y' must be numeric")
  if (length(x) != length(y))
    stop("'x' and 'y' must have the same number of observations")
  if (!is.numeric(lowerk))
    stop("argument 'lowerk' must be numeric")
  if (lowerk<=0)
    stop("argument 'lowerk' must be greater than 0")
  if (!is.numeric(lowerknu))
    stop("argument 'lowerknu' must be numeric")
  if (lowerknu<=0)
    stop("argument 'lowerknu' must be greater than 0")
  if (!is.numeric(upperk))
    stop("argument 'upperk' must be numeric")
  if (upperk<=0)
    stop("argument 'upperk' must be greater than 0")
  if (!is.numeric(upperknu))
    stop("argument 'upperknu' must be numeric")
  if (upperknu<=0)
    stop("argument 'upperknu' must be greater than 0")
  
  kappa_cv<-optimize(function(kappa)loglik_cv(x,y,kappa,startv),lower=lowerk,upper=upperk,maximum=TRUE)$maximum
  nucv_seq<-seq(lowerknu,upperknu,length=50)
  n_cv<-sapply(1:length(nucv_seq),function(i){my_fun_loglik_cv_nu(x,y,startv,startv2,kappa_cv,nucv_seq[i])})
  nu_cv<-nucv_seq[which.max(n_cv)]	
 
  sm_param<-c(kappa_cv,nu_cv)
  
  return(sm_param)
  
}

