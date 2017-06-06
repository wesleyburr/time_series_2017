#
#  Implementation of ACF (Autocorrelation Function)
#
myacf <- function(x, maxLag = 20) {

  N <- length(x)
  stopifnot(maxLag > 0)
  stopifnot(is.numeric(x), is.numeric(maxLag), is.vector(x), length(maxLag) == 1, length(which(is.na(x))) == 0,
            length(x) > 1, length(x) > maxLag)

  tol <- 1e-8 # tolerance level to be used 
  if (min(abs(c(maxLag %% 1, maxLag %% 1 - 1))) >= tol) {
    maxLag <- floor(maxLag)
    warning("Your maxLag is not an integer within 8 decimal places; the floor of the maxLag you entered has been used instead.")
  }

  # Computation
  xbar <- mean(x)
  xt_xbar <- x - xbar
  denom <- (xt_xbar)^2
  r_k <- vector(length = maxLag + 1)

  r_k[1] <- 1
  for (k in 1:maxLag)
  {
    first_part <- xt_xbar[1:(N-k)]
    second_part <- xt_xbar[(k+1):N]
    r_k[k + 1] <- as.numeric((t(first_part) %*% second_part)) / sum(denom)
  }
  names(r_k) <- paste("Lag", 0:maxLag, sep = "")

  r_k
}

