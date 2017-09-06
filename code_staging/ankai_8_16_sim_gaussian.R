#' @title the Innovation Algoritghm and simulation of Gaussian process
#'
#' @description The function ina() computes coefficient of the prediction of n+1 step as linear sun of 0 to n steps.
#' @description The input require autocovariance matrix
#' @description As return, the theta COLUMN index i representing the prediction of x_i (theta_{i-1,-}) and v is a error for corresponding prediction from v_0 to v_{n-1}
#'   
#' 
#'
#'
#' @param k autocovariance matrix for given time series
#'
#' @return v       (error)
#' @return theta     (prediction coefficients)
#'
#' @author Ankai Liu, \email{ankailiu@trentu.ca}
#' @references Brockwell, P. J., & Davis, R. A. (2013). Time series: theory and methods. New York, NY: Springer. p172 proposition 5.2.2
#' @seealso None
#' @keywords Inovation algorithm, linear prediction
#'
#' @examples
#' 
#' #example from book
#' 
#' n <- matrix(c(0.125, 0.05,  0,     0,     0,
#'               0.05,  0.125, 0.05,  0,     0,
#'               0,     0.05,  0.125, 0.05,  0,
#'               0,     0,     0.05,  0.125, 0.05,
#'               0,     0,     0,     0.05,  0.125),
#'               nrow = 5)
#' result <- ina(n)
#' result
#' 
#' #the actual solution
#'  
#' v <- (1+1/4)/100
#' for (i in 2:5) {
#'   v[i] <- (1+1/4-1/v[i-1]/4/100)/100
#' }
#' v-result$v
#' 
#' # another example
#' 
#' r <- 10
#' n <- matrix(rnorm(r*r),nrow = r)
#' ina(n)
#'




#
# since method #1 is not acceptable for computation. we thus use method #2
# the two inputs would be both strings with the coefficient.
# the output would be just a string with lag=c(0:20)
#


#
# The innovations algorithm
#
#
# k is the autocovariance matrix
# for given n by n matrix, we can only find prediction of X_n by \{X_1:X_{n-1}\}
#

#
# as return, the theta COLUMN index i representing the prediction of x_i (theta_{i-1,-}) and 
#     v is a vector from v_0 to v_{n-1}
#

ina <- function(k){
  stopifnot(is.numeric(k))
  stopifnot(length(k[1,]) == length(k[,1]))
  stopifnot(det(k) != 0)
  
  m <- length(k[1,]) 
  stopifnot(m > 1)
  
  
  v <- c(k[1,1])
  theta <- matrix(0, nrow = m, ncol = m)
  
  for (i in 1:(m - 1)) {  # theta_{i,*} and v_(i+1) in each step
    theta[i,i + 1] <- 1/v[1] * k[i + 1, 1] # j=0
    for (j in (i - 1):0) {
      theta[i - j,i + 1] <- 1/v[j + 1] * (k[i + 1, j + 1] - 
                                 sum((theta[,(j + 1)])[j:1] * (theta[,(i + 1)])[i:(i - j + 1)] * v[1:j])
                               ) 
    }
    
    v_t <- k[i + 1, i + 1] - sum((theta[,(i + 1)])[i:1]^2 * v[1:i])
    v <- c(v, v_t)   #v_i
  }
  
  return(
    list("v" = v, "theta" = theta)
  )
}

#1/10*(1+1/4)
#1/20
#n <- matrix(c(0.125, 0.05,  0,     0,     0,
#              0.05,  0.125, 0.05,  0,     0,
#              0,     0.05,  0.125, 0.05,  0,
#              0,     0,     0.05,  0.125, 0.05,
#              0,     0,     0,     0.05,  0.125),
#            nrow = 5)
#n
#result <- ina(n)
#result$v
#result
#
#v <- (1+1/4)/100
#for (i in 2:5) {
# v[i] <- (1+1/4-1/v[i-1]/4/100)/100
#}
#v
#result$theta
#v-result$v

r <- 10
n <- matrix(rnorm(r*r),nrow = r)
ina(n)

#
# simulation of a gaussian process
#

#
# the input is an autocovariance matrix k and 
#     np be the # of observation points we want to simulated
#


#' @title Simulation of a Gaussian Process
#'
#' @description The function sgp() computes Simulation of a Gaussian process for given autocovariance matrix and desired number of steps.
#' @description the input would be autocovariance matrix and the number of time steps.
#' @description the output the estimated gaussian process.
#'
#' @param k autocovatiance matrix
#' @param np the number of time steps we wnat to use
#'
#' @return coe       (coefficient of estimation from Z_1 to Z_n)
#'
#' @note this algorithm based on inovation algorithm and require the v (error) to be positive.
#'
#' @author Ankai Liu, \email{ankailiu@trentu.ca}
#' @references Brockwell, P. J., & Davis, R. A. (2013). Time series: theory and methods. New York, NY: Springer. p271 8.16
#' @seealso None
#' @keywords estimation, guassian process
#'
#' @examples
#'
#' n <- matrix(c(0.125, 0.05,  0,     0,     0,
#'               0.05,  0.125, 0.05,  0,     0,
#'               0,     0.05,  0.125, 0.05,  0,
#'               0,     0,     0.05,  0.125, 0.05,
#'               0,     0,     0,     0.05,  0.125),
#'             nrow = 5)
#' sgp(n,4)
#' 
#' r <- 10
#' k <- matrix(rnorm(r*r),nrow = r)
#' sgp(k, r-1)
#'
#' @export
#' @importFrom grDevices rgb2hsv
#' @importFrom graphics par plot rect text
#'

sgp <- function(k, np) {
  stopifnot(is.numeric(k), is.numeric(np))
  stopifnot(length(k[1,]) == length(k[,1]))
  stopifnot(det(k) != 0)
  stopifnot(np >= 2)
  stopifnot(length(k[1,]) >= np + 1)
  
  if(ceiling(np) != np) {
    lag <- ceiling(np)
  }
  theta <- ina(k)$theta[,np]
  v <- ina(k)$v
  
  if (sum(n <= 0) > 0) {
    print("v=")
    print(v)
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    warning("negative v, see inovation algorithm")
  }
  
  #print(theta[(np - 1):1])
  #print(v[1:(np)])
  coe <- theta[(np - 1):1] * v[1:(np - 1)]^(1/2)  
  
  #
  # *** here we might have negatice v and cause problems
  #
  
  coe <- c(coe, v[np]^(1/2))
  return(coe)
}



#r <- 10
#k <- matrix(rnorm(r*r),nrow = r)
#sgp(k, r-1)
#
#n <- matrix(c(0.125, 0.05,  0,     0,     0,
#              0.05,  0.125, 0.05,  0,     0,
#              0,     0.05,  0.125, 0.05,  0,
#              0,     0,     0.05,  0.125, 0.05,
#              0,     0,     0,     0.05,  0.125),
#           nrow = 5)
#sgp(n,4)
