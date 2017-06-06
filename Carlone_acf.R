#
# Carlone Scott
# Time Series Assigment
# Chapter 2
# Prof Wesley Burr
#

# ACF Function

#
#  Making an ACF Function
#
manual_acf <- function(x, maxLag = 20) {
  # Sanity Checks?!
  stopifnot(is.numeric(x), is.numeric(maxLag), is.vector(x), length(x) > 1)
  stopifnot(length(which(is.na(x))) == 0)
  
  # Vector shorter than maxLag?
  
  # Cast the maxLag to integer ... then see if it changed ... if so, warning()
  
  If (maxLag - round(maxLag) !=0){
    maxLag <- round(maxLag)
  }
  stopifnot(length(x) > abs(maxLag))
  
  warning("OMG YOU SCREWED UP") 
  stop("OMG YOU DID EVEN WORSE")
  
  # Compute ACF
  
  xbar <- sum(x)/n
  rho <- (sum((x[-c(1:maxLag)] - xbar)*(x[-c(n:(n-maxLag+1))] - xbar)))/(sum((x - xbar)^2))
  
  
  # Return
  acf
}

# Calls
acf(x)
acf(x, whichForm = "version2")

# acf(x)$acf == manual_acf

