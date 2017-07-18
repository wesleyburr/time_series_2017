#' @title the Innovation Algoritghm and simulation of Gaussian process
#'
#' @description 
#' 
#' 
#' @param 
#'
#' @return 
#'
#' @note This function STRONGLY REQUIRE arma process to be causal.
#'
#' @author Ankai Liu, \email{ankailiu@trentu.ca}
#' @references Brockwell, P. J., & Davis, R. A. (2013). Time series: theory and methods. New York, NY: Springer.
#' @seealso None
#' @keywords 
#' 
#' @examples
#' 
#' @export
#' @importFrom grDevices rgb2hsv
#' @importFrom graphics par plot rect text
#'


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

1/10*(1+1/4)
1/20
n <- matrix(c(0.125, 0.05,  0,     0,     0,
              0.05,  0.125, 0.05,  0,     0,
              0,     0.05,  0.125, 0.05,  0,
              0,     0,     0.05,  0.125, 0.05,
              0,     0,     0,     0.05,  0.125),
            nrow = 5)
n
result <- ina(n)
result$v

v <- (1+1/4)/100
for (i in 2:5) {
  v[i] <- (1+1/4-1/v[i-1]/4/100)/100
}
v
result$theta
v-result$v

r <- 10
n <- matrix(rnorm(r*r),nrow = r)
ina(n)$v

#
# simulation of a gaussian process
#

#
# the input is an autocovariance matrix k and 
#     np be the # of observation points we want to simulated
#

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
  
  
  print(theta[(np - 1):1])
  #print(v[1:(np)])
  coe <- theta[(np - 1):1] * v[1:(np - 1)]^(1/2)  
  
  #
  # *** here we might have negatice v and cause problems
  #
  
  coe <- c(coe, v[np]^(1/2))
  return(coe)
}



r <- 10
k <- matrix(rnorm(r*r),nrow = r)
sgp(k, r-1)

n <- matrix(c(0.125, 0.05,  0,     0,     0,
              0.05,  0.125, 0.05,  0,     0,
              0,     0.05,  0.125, 0.05,  0,
              0,     0,     0.05,  0.125, 0.05,
              0,     0,     0,     0.05,  0.125),
            nrow = 5)
sgp(n,4)
