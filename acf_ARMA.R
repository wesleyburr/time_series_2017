#acvf for ARMA process

#
# since method #1 is not acceptable for computation. we thus use method #2
# the two inputs would be both strings with the coefficient.  
# the output would be just a string with lag=c(0:20) 
#


acvf_arma <- function( phi, theta, lag = 20 , output) {
  stopifnot(is.numeric( phi ), is.numeric( theta ), 
            is.vector( phi ), is.vector( theta ), 
            length( phi ) > 0, length( theta ) > 0, 
            is.character( output ))
  
  if (phi[1] == 0) {
    stop( "leading coefficient of PHI is 0")
  }
  if (theta[1] == 0) {
    stop( "leadning coefficient of THETA is 0")
  }
  
  phi <- phi / phi[1]
  theta <- theta / theta[1]
  #
  # we first find psi=theta/phi
  #
  p <- length( phi )
  q <- length( theta )
  
  phi_1 <- c( phi, rep(0,40) )
  theta_1 <- c( theta, rep(0,40) )
  psi <- rep(0,40)
  psi[1] <- theta_1[1]
  for (i in 2:40) {
    psi[i]=theta_1[i] - sum(psi[1:(i - 1)] * phi_1[i:2])
  }
  #
  # method 2 find lag 1:p 
  # 
  
  #
  # Left hand side
  #
  phi_2 <- c(rep(0,p),phi,rep(0,2 * p))
  l <- matrix(c(phi_2[(p + 1):(2 * p)],rep(0,(p)^2 - p)),nrow=p,ncol=p,byrow=FALSE)
  for (i in 1:p) {
    for (j in 2:p) {
      l[i,j] <- phi_2[p + i - (j - 1)]+phi_2[p + i + (j - 1)]
      #print(phi_2[p - 1 + i - (j - 1)])
      #print(phi_2[p - 1 + i + (j - 1)])
    }
  }
  r <- rep(0,p)
  
  #
  # right hand side
  # 
  
  for (i in 1:21) {
    if ( i <= p) {
      if (i <= q) {
        r[i] <- sum(psi[1:(q - i + 1)] * theta[i:q])
        }
    }
    else {
      r[i] <- 0
    }
  }
  gam <- solve(l,r[1:p])
  
  #
  # method 3, we use existence gamma to find the rest
  #
  if ( p < 21 ) {
    for (i in (p + 1):21) {
      g <- -sum(phi[-1] * gam[(i - 1):(i - p + 1)])
      gam <- c(gam, g)
    }
  }
  #
  # output
  #
  rho <- gam / gam[1]
  output <- data.frame(lag=0:20,gamma=gam, rho=rho)
#   if (output == "acvf") {
#     return (gam)
#   }
#   else {
#     return (rho) 
#   }
  return(output)
}


#
# test acf
#
phi <- c(1, -1, 1 / 4)
theta <- c(1, 1)

testgamma <- rep(0, 21)
for (i in 1:21) {
  testgamma[i] <- (32 / 3 + 8 * (i - 1)) * 2^(- i + 1)
}
testrho <- testgamma / testgamma[1]
test1 <- data.frame(lag=0:20, gamma=testgamma,rho=testrho)
test1

test2 <- acvf_arma(phi, theta, 20, "acf")
library('ggplot2')
ggplot(test2, aes(x=lag,y=gamma))+geom_bar(stat = "identity")+geom_bar(data=test1, aes(x=lag,y=gamma),color="red",stat="identity")
ggplot(test2, aes(x=lag,y=rho))+geom_bar(stat = "identity")+geom_bar(data=test1, aes(x=lag,y=rho),color="red",stat="identity")
test2

test1-test2

#
# test acvf
#
testgamma
barplot(testgamma)

test2 <- acvf_arma(phi, theta, 20, "acvf")
barplot(test2)
test2

barplot(testgamma-test2)


#
# test psy
#

testpsy <- rep(0,20)
for (i in 1:20) {
  testpsy[i] <- (1+3*(i-1))*2^(-i+1)
}
testpsy-psi

#
# else
#
phi <- c(0, 1, -1, 1 / 4)
theta <- c(1, 1, 0, 0, 0)
acvf_arma(phi, theta, 20, "acvf")

phi <- c( 1, -1, 1 / 4)
theta <- c( 0, 1, 1, 0, 0, 0)
acvf_arma(phi, theta, 20, "acvf")
