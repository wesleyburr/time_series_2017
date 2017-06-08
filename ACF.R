#
# myacf
# mw

my_acf <- function(x, maxlag) 
{ #Start
  
  #Stop(1)------------------------------------------------------------------
  stopifnot(is.numeric(x), is.numeric(maxlag), is.vector(x), length(x)>1)
  stopifnot(length(which(is.na(x))) == 0)
  #-------------------------------------------------------------------------
  # Cast the maxlag to integer .. then see if it changed .. if so, warning()
  
  if (maxlag-round(maxlag)!=0){
    warning('wrong input type')
    maxlag  <- round(maxlag)
  }
  
  #Stop(2)------------------------------------------------------------------   
  stopifnot(length(x)>maxlag)
  if (maxlag == 0){
    return(1)
  }
  #-------------------------------------------------------------------------
  
  if (maxlag == 0){
    return(1)
  }
  
  n <- length(x)
  m <- mean(x)
  
  for(i in 1:maxlag){
    rho = (sum((x[-c(1:i)]-m)*(x[-c(n:(n-i+1))]-m)))/(sum((x-m)*(x-m)))
    print(rho)}
  
}


#Testing
# call
y <- c(1:5)
y <- rnorm(1000)
my_acf(y,2)
acf(y,2)
