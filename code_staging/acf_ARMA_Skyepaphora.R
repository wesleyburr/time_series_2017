ARMA_acf <- function(phi,theta){
  stopifnot(is.numeric( phi ), is.numeric( theta ), 
            is.vector( phi ), is.vector( theta ), 
            length( phi ) > 0, length( theta ) > 0)
  ##orders (should I be concerned about the way R uses indeces? The textbook uses indeces starting at 0 and it throws me off but I tried to work arund it):
  phi <- -(phi)
  p <- length(phi)
  q <- length(theta)
  lo_ord <- min(p,q)
  ##now look for psi values
  psi <- c(theta[1],rep(0,(lo_ord-1)))
  for (i in (2:lo_ord)){
    phi2 <- rep(phi[2:i])
    psi2 <- rep(psi[(i-1):1])
    psi[i] <- theta[i] + t(phi2)%*%psi2
  }
  #create a matrix of phi's, a vector for the right-hand side, and solve
  phivec <- rep(0,p)
  bigphi <- rep(0,(p^2))
  phi_w_0 <- c(rep(0,p),phi,rep(0,2*p))
  RS <- rep(0,p)
  RS
  psi3 <- rep(0:p)
  for (i in 1:p) {
    for (j in 2:p){ 
      phivec[1] <- phi[i]
      phivec[j] <- phi_w_0[p + i - (j - 1)]+phi_w_0[p + i + (j - 1)]
    }   
    bigphi[(((i-1)*p)+1):(i*p)] <- phivec
    theta2 <- c(theta[i:q])
    if (i <= q){
      psi3 <- c(psi[1:(q-i+1)])
      RS[i] <- t(theta2)%*%psi3
    }
    else{
      RS[i] <- 0
    }
  }
  phi_mat <- matrix(bigphi,nrow = p, ncol = p, byrow = TRUE)
  gamma <- solve(phi_mat,RS)
  return(-(gamma))
}
