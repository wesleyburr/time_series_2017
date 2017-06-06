#
# ACF function (sample autocovariance function)
# Ankai Liu	

#
# autocovariance
#
myacv <- function(x, maxlag) {
  stopifnot(is.numeric(x), is.numeric(maxlag), is.vector(x), length(x)>1)
  stopifnot(length(which(is.na(x))) == 0)
  
  # Cast the maxlag to integer .. then see if it changed .. if so, warning()
  if (maxlag-round(maxlag)!=0) {
    warning('wrong input type')
    maxlag  <- round(maxlag)
  }
  stopifnot(length(x)>abs(maxlag))
  
  maxlag <- abs(maxlag)
  n <- length(x)
  mean <- sum(x)/n
  #
  # for maxlag=0 we return variance
  #
  if (maxlag==0) {
    return((sum((x-mean)*(x-mean))/n))
  }
  gamma <- sum((x[-c(1:maxlag)]-mean)*(x[-c(n:(n-maxlag+1))]-mean))/n
  return(gamma)
}

# call
x <- c(1:5)
x <- rnorm(1000)
myacv(x,0)
myacv(x,2)
myacv(x,1)
(built_in_cv <- acf(x, type = "covariance"))
#
# autocorelation
#
myacf <- function(x, maxlag) {
  stopifnot(is.numeric(x), is.numeric(maxlag), is.vector(x), length(x)>1)
  stopifnot(length(which(is.na(x))) == 0)
  
  # Cast the maxlag to integer .. then see if it changed .. if so, warning()
  if (maxlag-round(maxlag)!=0) {
    warning('wrong input type')
    maxlag  <- round(maxlag)
  }
  stopifnot(length(x)>abs(maxlag))
  
  #
  # because for acf lag 0 always return 1
  #
  if (maxlag == 0) {
    return(1)
  }
  maxlag <- abs(maxlag)
  n <- length(x)
  mean <- sum(x)/n
  rho <- (sum((x[-c(1:maxlag)]-mean)*(x[-c(n:(n-maxlag+1))]-mean))/n)/(sum((x-mean)*(x-mean))/n)
  return(rho)
}
# call
x <- c(1:5)
x <- rnorm(1000)
myacf(x,0)
myacf(x,1)
myacf(x,2)
myacf(x,3)
myacf(x,-3.8)
(built_in_cr <- acf(x, type = "correlation") )




