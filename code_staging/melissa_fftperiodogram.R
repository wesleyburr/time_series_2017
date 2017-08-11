fftperiodogram <- function(x, seconds) {
  # Catches and warnings 
  stopifnot(is.numeric(x), length(x) > 0, seconds >= 0, is.numeric(seconds))
  warning(if(length(x) %% 2 == 1) {
    print("Must have even number of points; first point was removed")
  })
  
  # Ensure length of time series is an even number; requirement for R's fft function
  n <- length(x)
  if(n %% 2 == 1) {
    x <- x[0:(n-1)]
    n <- length(x)
  }
  
  # Calculate the discrete fourier transform
  dft <- abs(fft(x)/sqrt(n))^2
  w <- n/2 + 1
  P <- (4/n)*dft[1:w]
  freq <- (0:(n/2))/n
  
  # Calculate the frequency at which the maximum value occurs in the periodogram and print out the corresponding period
  maximumfreq <- max(P[2:length(P)])
  location <- which(P[1:length(P)] >= maximumfreq)
  period <- 1/freq[location[2]]
  sentence <- paste("There appears to be a pattern around every", period*seconds, "seconds, which is", period*seconds/(86400*365), "years", sep = " ")
  
  # Finally, plot the periodogram
  periodogram_plot <- plot(freq, P, type = "h", lwd = 2, xlab = "frequency", ylab = "I(lambda)", xlim = c(0.019, 0.5), ylim = c(0, maximumfreq + 10), 
                           main = "Periodogram")
  
  # Make a list of the things that the function should return, and return it  
  periodogram_return <- list(sentence, periodogram_plot)
  return(periodogram_return)
}


