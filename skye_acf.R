#This is {x_t}_{t\in\mathcal{T}} where \mathcal{T} is {1,...,n} and yes I'm using little n's
sample_acf <- function(x,k = 1){
  #Catch bad things
  stopifnot(is.numeric(x), is.numeric(k), is.vector(x), length(x) > 1)
  stopifnot(length(which(is.na(x))) == 0)
  stopifnot(k < length(x))
  #This is the sample mean \overbar{x}
  xbar <- mean(x)
  #This is (x_t - \overbar{x})
  y <- x - xbar
  #This is the sample variance \hat{sigma}^2
  var_hat <- sum(y^2)
  #The sum of y times w 
  dot <- t(y[1:(n-k)]) %*% y[(k+1):n]
  #Standardized (divide by sample variance)
  rho_k <- dot/var_hat
  return(rho_k)
}
n <- 10
x = rnorm(n)
sample_acf(x,1)
acf(x)$acf[2]
