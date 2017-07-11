#' Durbin Levinson Calculator
#'
#' Takes a vector of autocovariances, gamma, and a vector of realizations, X, and computes a prediction for the next realization, using the Durbin Levinson algorithm.
#'
#' @author Melissa Van Bussel <melissavanbussel@trentu.ca>
#'
#' @keywords Durbin Levinson, Durbin, Levinson, ARMA
#'
#' @param gamma A numeric vector of autocovariances, beginning at gamma(0) and ending at gamma(n)
#' @param X A numeric vector of realizations, beginning at X_1 and ending at X_n
#'
#' @return The prediction for X_(n+1)
#'
#' @examples
#' # vectors
#' X <- c(6, 7, 8, 9) # length of X should be one less than length of gamma
#' gamma <- c(5, 4.3, 3.2, 2.7, 1.6)
#'
#' # Use function
#' durbin_levinson(gamma = gamma, X = X)
#'
#' @references Brockwell, P. J., & Davis, R. A. (2013). Time series: Theory and Methods. New York, NY: Springer.
#'
#' @export

durbin_levinson <- function(gamma, X) {
  # Catches
  stopifnot(is.numeric(X), is.numeric(gamma), length(X) == (length(gamma)-1))

  # Calculations
  n <- length(X)
  rho <- gamma/gamma[1]
  rho <- rho[-1]

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
  X_next <- sum(phi[,n] * rev(X))
  return(X_next)
}
