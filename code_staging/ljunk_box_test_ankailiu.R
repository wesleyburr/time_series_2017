
#' @title Ljung-Box test for AR proces
#'
#' @description for given observation of AR process and estimation of coefficients
#' @description the Ljung-Box test test the lack of fit for our estimation
#'
#' @param phi estimated coefficient of AR process
#' @param X observation of AR process
#' @param m number if lag of autocovariance of residual we want to use
#'
#' @return Q        (Q)
#' @return cut off  (chi-square with corresponding degree of freedom)
#'
#' @note if our estimation is $X_t-\sum_{j=1}^p \phi_j X_{t-j}=Z_t$, $\phi=\{\phi_1, \phi_2,\ cdots, \phi_p\}$ 
#'
#' @note if $Q>cut off$, we reject that we get the true coefficients. 
#'       if $Q<cut off$, we cannt reject that we get the true coefficients 
#'
#' @author Ankai Liu, \email{ankailiu@trentu.ca}
#' @references Box, G. E., $\&$ Pierce, D. A. (1970). Distribution of Residual Autocorrelations in Autoregressive-Integrated Moving Average Time Series Models. Journal of the American Statistical Association, 65(332), 1509. doi:10.2307/2284333
#' @references Ljung, G. M., $\&$ Box, G. E. (1978). On a Measure of Lack of Fit in Time Series Models. Biometrika, 65(2), 297. doi:10.2307/2335207#' @seealso None
#' 
#' @keywords Ljung-Box test, lack of fit
#'
#' @examples
#' # Define time series data in a vector, x
#' 
#' n <- 1000
#' x <- arima.sim(n=n, list(ar = c(1,-1,1/4)), sd = sqrt(2))
#' acv <- c(myacv(x,0),myacv(x,1),myacv(x,2),myacv(x,3),myacv(x,4),myacv(x,5),myacv(x,6))
#' plot(abs(acv), type="l")
#' 
#' phi1 <- c(1/2,-1,1)
#' ljungbox(x,phi1,10)
#' ljungbox(x,phi1,20)
#' ljungbox(x,phi1,30)
#' 
#' phi2 <- c(1/4,-1)
#' ljungbox(x,phi2,10)
#' ljungbox(x,phi2,20)
#' ljungbox(x,phi2,30)
#' @export
#' @importFrom grDevices rgb2hsv
#' @importFrom graphics par plot rect text
#'

ljungbox <- function(x, phi, m) {
  
  # we want n way larget than p and m
  
  p <- length(phi)
  n <- length(x)
  a <- c()
  r <- c()
  
  for (i in (p + 1):n) {
    a[i - p] <- x[i] - sum(phi * x[(i - p):(i - 1)])
  }
  k <- length(a)
  for (i in 1:(k - 1)) {
    r[i] <- sum(a[(i + 1):k] * a[1:(k - i)])/sum(a * a)
  }
  
  Q <- n * (n + 2) * sum(r[1:m]^2/c((n - 1):(n - m)))
  cutoff <- qchisq(0.95, df = (m - p))
  return(data.frame("Q"= Q , "cut off" = cutoff))
}



#########test
#n <- 1000
#x <- arima.sim(n=n, list(ar = c(1,-1,1/4)), sd = sqrt(2))
#acv <- c(myacv(x,0),myacv(x,1),myacv(x,2),myacv(x,3),myacv(x,4),myacv(x,5),myacv(x,6))
#plot(abs(acv), type="l")
# 
#x
#
#phi1 <- c(1/2,-1,1)
#ljungbox(x,phi1,10)
#ljungbox(x,phi1,20)
#ljungbox(x,phi1,30)
#
#phi2 <- c(1/4,-1)
#ljungbox(x,phi2,10)
#ljungbox(x,phi2,20)
#ljungbox(x,phi2,30)

