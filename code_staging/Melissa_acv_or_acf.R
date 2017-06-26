acv_or_acf <- function(x, form = "acf") {
  q <- length(x)
  # Catches
  stopifnot(is.numeric(x), form %in% c("acf", "acv"), q <= 20)
  
  # Computations
  acf_vector <- vector(length = q)
  acv_vector <- vector(length = q)
  
  for (i in 1:q) {
    acv_vector[i] <- sum(x[1:(q-i+1)]*x[i:q])
  }
  
  
  # What to return
  if (form == "acf")
  {
    acf_vector <- acv_vector/acv_vector[1]
    return(acf_vector)
  }
  else {
    return(acv_vector)
  }
}