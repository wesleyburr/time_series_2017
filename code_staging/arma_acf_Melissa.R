#' ACF calculator for ARMA Processes
#'
#' Takes a phi vector and theta vector and returns the acv or acf of the ARMA process.
#'
#' @author Melissa Van Bussel <melissavanbussel@trentu.ca>
#'
#' @keywords ARMA, ACF, AR, MA, ACV, Autocovariance, Autocorrelation
#'
#' @param phi A numeric vector of the phi coefficients.
#' @param theta A numeric vector of the theta coefficients.
#' @param which_form Optional argument: either acv or acf; indicates which function is returned to the user.
#'
#' @return Either the acv or the acf of the ARMA process.
#'
#' @examples
#' # Create phi and theta vectors
#' phi <- c(1, -1.3, 0.4)
#' theta <- c(1, 1)
#'
#' # Use function
#' arma_acf(phi, theta, which_form = "acv")
#'
#' @references Brockwell, P. J., & Davis, R. A. (2013). Time series: Theory and Methods. New York, NY: Springer.
#'
#' @export


arma_acf <- function(phi, theta, which_form = "acf") {
  # Catches
  stopifnot(is.numeric(phi), is.numeric(theta), length(theta) <= 20, length(phi) <= 20)

  # Computations
  # How many lags to calculate
  p <- length(phi)
  q <- length(theta)
  max_value <- max(p - 1, q)
  # Make all 3 vectors length 20
  phi <- c(phi, rep(0, 20-p))
  theta <- c(theta, rep(0, 20-q))
  psi <- rep(0, 20)
  for (i in 1:max_value) {
    psi[i] <- (theta[i] - sum(psi*rev(phi[0:max_value]))) / phi[1]
  }
  psi <- psi[0:max_value]

  # Solve matrix to get acv
  rhs <- vector(length = max_value + 1)
  for (j in 1:max_value) {
    rhs[j] <- sum(psi[1:(max_value + 1 - j)])
  }

  # Remove extra zeroes now that we are done with them
  phi <- phi[1:p]

  # LHS (start with filling 1st row and 1st column)
  lhs <- matrix(0:0, ncol = p, nrow = p)
  lhs[1,] <- phi
  lhs[,1] <- phi

  # Matrix of repeated rows of phi
  phi_mat <- matrix(ncol = p, nrow = p)
  for (i in 1:p) {
    phi_mat[i,] <- phi
  }

  # Matrix with only the diagonals you care about
  for (i in 1:(p-2)) {
    diag_mat <- matrix(0:0, ncol = p, nrow = p)
    diag(diag_mat[,-(1:i)]) <- diag(phi_mat[, -(1:i)])
    diag(diag_mat[-(1:i),]) <- diag(phi_mat[-(1:i),])
    diag_mat_vector <- apply(diag_mat, 1, sum) # 1 = dimcode, corresponds to row
    lhs[,(i+1)] <- diag_mat_vector
  }
  # can't take diagonal if dim is less than 2, so bottom right corner is done manually:
  lhs[p,p] <- phi[1]

  acv_1 <- solve(lhs, rhs)
  acf_1 <- acv_1/acv_1[1]

  if(which_form == "acv") {
    names(acv_1) <- paste("Lag", 0:(length(acv_1) - 1), sep = "")
    return(acv_1)
  } else {
    names(acf_1) <- paste("Lag", 0:(length(acv_1) - 1), sep = "")
    return(acf_1)
  }
}
