x <- seq(1:30)
acf_k <- function(x, maxLag = 20)
{
  N <- length(x)
  
  # Catches
  stopifnot(is.numeric(x), is.numeric(maxLag), is.vector(x), length(maxLag) == 1, length(which(is.na(x))) == 0)
  
  if (!(N > 1)) {
    stop("Your time series, x, must be at least two points.")
  } 
  
  if (is.vector(maxLag) & length(maxLag) > 1) {
    stop("maxLag must be a numeric (integer); should not be a vector.")
  }
  
  # Warnings
  if (maxLag > N) {
    maxLag <- N-1
    warning("Your maxLag is greater than the length of your time series; it has been truncated to the length of the time series minus one.")
  }
  
  tol <- 1e-8 # tolerance level to be used 
  if (min(abs(c(maxLag%%1, maxLag%%1-1))) >= tol) {
    maxLag <- floor(maxLag)
    warning("Your maxLag is not an integer within 8 decimal places; the floor of the maxLag you entered has been used instead.")
  }
  
  # Computation
  xbar <- mean(x)
  xt_xbar <- x - xbar
  denom <- (xt_xbar)^2
  r_k <- vector(length = maxLag)
  
  for (k in 1:maxLag)
  {
    first_part <- xt_xbar[1:(N-k)]
    second_part <- xt_xbar[(k+1):N]
    r_k[k] <- as.numeric((t(first_part)%*%second_part))/sum(sigma2) 
  }
  r_k <- c(1, r_k)
  names(r_k) <- paste("Lag", 0:maxLag, sep = "")
  return(r_k)
}