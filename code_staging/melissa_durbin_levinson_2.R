#' Durbin Levinson Calculator Version 2
#'
#' Takes a vector of realizations, x, and uses it to create a vector of sample autocorrelations. The sample autocorrelations are then used to predict the next realization, via the Durbin-Levinson algorithm.
#'
#' @note If the length of x is N, the length of the sample autocorrelation vector will be N-1. This means that the Durbin-Levinson algorithm is calculating a prediction for the Nth value in the time series, not the (N+1)th.
#'
#' @author Melissa Van Bussel <melissavanbussel@trentu.ca>
#'
#' @keywords Durbin Levinson, Durbin, Levinson, ARMA
#'
#' @param x A numeric vector of realizations, beginning at x_1 and ending at x_N
#'
#' @return The prediction for x_N
#'
#' @examples
#' # vector
#' x_vector <- c(6, 7, 8, 9)
#'
#' # Use function
#' durbin_levinson_2(x = x_vector)
#'
#' @references Brockwell, P. J., & Davis, R. A. (2013). Time series: Theory and Methods. New York, NY: Springer.
#'
#' @export

durbin_levinson_2 <- function(x) {

  N <- length(x)

  # Catches
  stopifnot(is.numeric(x), is.vector(x), length(which(is.na(x))) == 0)

  if (!(N > 1)) {
    stop("Your time series, x, must be at least two points.")
  }

    maxLag <- N-1

  # Computation of sample acv
  xbar <- mean(x)
  xt_xbar <- x - xbar
  denom <- (xt_xbar)^2
  rho <- vector(length = maxLag)

  for (k in 1:maxLag)
  {
    first_part <- xt_xbar[1:(N-k)]
    second_part <- xt_xbar[(k+1):N]
    rho[k] <- as.numeric((t(first_part)%*%second_part))/sum(denom)
  }

# Computation of phi matrix from Durbin-Levinson algorithm
n <- length(rho)
phi <- matrix(0:0, ncol = n, nrow = n)
phi[1,1] <- rho[1]

for (m in 2:n) {
  for (j in 2:n) {
    if (m == j)
    {
      phi[m, m] <- (rho[m] - phi[1:(m-1), (m-1)] %*% rev(rho[-(m:n)])) / (1 - phi[1:(m-1),(m-1)] %*% rho[-(m:n)])
    }
    else
    {
      i <- 1:(j-1)
      k <- rev(i)
      for (i in i) {
        phi[i, j] <- phi[i, (j-1)] - phi[j,j]*phi[k[i], (j-1)]
      }
    }
  }
}
# for some reason it works exactly when you run this block twice... not sure what's up with that...
for (m in 2:n) {
  for (j in 2:n) {
    if (m == j)
    {
      phi[m, m] <- (rho[m] - phi[1:(m-1), (m-1)] %*% rev(rho[-(m:n)])) / (1 - phi[1:(m-1),(m-1)] %*% rho[-(m:n)])
    }
    else
    {
      i <- 1:(j-1)
      k <- rev(i)
      for (i in i) {
        phi[i, j] <- phi[i, (j-1)] - phi[j,j]*phi[k[i], (j-1)]
      }
    }
  }
}
x_next <- sum(phi[,n]*rev(x[2:(n+1)]))
return(x_next)
}
