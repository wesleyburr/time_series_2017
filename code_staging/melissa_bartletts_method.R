bartletts_method <- function(x, segments) {
  library(data.table)
  # Catches
  stopifnot(is.numeric(x), length(x) > 1)
  if (!is.ts(x)) {
    x <- as.ts(x)
    warning("x is not a time-series, conversion 'as.ts(x)' was applied.")
  }
  n <- length(x)
  # check if the time series can be divided evenly into the number of segments given, fix it if not.
  if (n %% segments != 0) {
    x <- c(x, rep(0, n %% segments))
    warning("length of time series not a multiple of segments, zeroes were added at end")
  }
  
  # Reassign length and find frequency
  name <- deparse(substitute(x))
  xfreq <- frequency(x)
  
  f <- 1:segments
  splitup <- split(x, f) # splits the time series into the number of specified segments 
  N <- length(splitup$`1`) # what is the length of each segment?
  w <- N/2 + 1 # time to calculate the periodograms for each of the segments 
  freq <- seq.int(from = xfreq/N, by = xfreq/N, length.out = w) # frequency for x axis 
  P <- list()
  meanvector <- vector(length = w)
  for(k in 1:w) {  
    for (i in 1:segments) {
      dft <- abs(fft(splitup[[i]])/sqrt(N))^2
      P[[i]] <- dft[1:w]
      Q <- transpose(P)
      meanvector[k] <- mean(Q[[k]])
    } 
  }
  
  plot(freq, log(meanvector), type = "l",
       xlab = "frequency", ylab = "log(averaged periodogram value)",
       main = paste("Averaged periodogram of", name,"dataset using Bartlett's Method", sep = " "),
       cex.main = 0.8)
}
