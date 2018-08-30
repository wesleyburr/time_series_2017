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
  gamma <- c()
  #
  # for maxlag=0 we return variance
  #
  for (i in c(0:maxlag)){
    if (i==0) {
      gamma <- (sum((x-mean)*(x-mean))/n)
    }
    gamma <- c(gamma,sum((x[-c(1:i)]-mean)*(x[-c(n:(n-i+1))]-mean))/n)
  }
  return(gamma)
}

# call
x <- rnorm(1000)
gamma <- myacv(x,0)
gamma
gamma <- myacv(x,2)
gamma
gamma <- myacv(x,10)
gamma


k <- c()
for (i in c(0:10)){
  k <- c(k,i)
}
k
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
  rho <- c()
  for (i in c(0:maxlag)){
    rho <- c(rho,(sum((x[-c(1:i)]-mean)*(x[-c(n:(n-i+1))]-mean))/n)/(sum((x-mean)*(x-mean))/n))
  }
  return(rho)
}
# call
x <- c(1:5)
x <- rnorm(1000)
rho <- myacf(x,0)
rho
rho <- myacf(x,1)
rho
rho <- myacf(x,2)
rho
rho <- myacf(x,30)
rho
rho <- myacf(x,-3.8)
rho
plot(rho,type = "h")



