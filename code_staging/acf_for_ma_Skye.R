ma_acv <- function(x, lag = 1){
  #fill vector with lag+1 many zeroes.
  y <- rep(0, lag + 1)
  if(lag > 0) {
    k <- min(length(x) - 1, lag)
    n <- length(x)
    for(i in 0:k) {
      y[(i+1)] <- t(x[1:(n-i)]) %*% x[(i+1):n]
    }
  } else {
    y[1] <- t(x) %*% x  
  }
  y
}

ma_acf <- function(x, lag = 0){
  z <- (ma_acv(x,lag)/ma_acv(x,0))
  return(z)
}