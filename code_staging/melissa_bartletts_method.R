bartletts_method <- function(x, segments) {
  # Catches
  stopifnot(is.numeric(x), length(x) > 1)
  if (!is.ts(x)) {
    x <- as.ts(x)
    warning("x is not a time-series, conversion 'as.ts(x)' was applied.")
  }
  n <- length(x)
  # check if the time series can be divided evenly into the number of segments given, fix it if not.
  if (n %% segments != 0) {
    x <- x[1:((floor(n / segments))*segments)]
    warning("length of time series not a multiple of segments, points were cut off at end")
  }
  f <- 1:segments
  splitup <- split(x, f) # splits the time series into the number of specified segments 
  N <- length(splitup$`1`) # what is the length of each segment?
  w <- N/2 + 1 # time to calculate the periodograms for each of the segments 
  freq <- (0:(N/2))/N
  P <- list()
  meanP <- matrix(data = NA, nrow = w, ncol = segments) 
  meanvector <- vector(length = w)
  for (k in 1:w) {
    for (i in 1:segments) {
      dft <- abs(fft(splitup[[i]])/sqrt(N))^2
      P[[i]] <- (4/N)*dft[1:w]
      meanP[k,i] <- P[[i]][k]
      meanvector[k] <- mean(meanP[k,])
    } }
  # calculate maximum frequency to choose ylim for the combined periodogram plot 
  maximumfreq <- max(meanvector[2:length(meanvector)])
  # finally, plot the average of all the periodograms
  combined_periodogram <- plot(freq, meanvector, ylim = c(0, maximumfreq*1.5), type = "h", xlim = c(0.02, 0.5), ylab = "average of I(lambda)", lwd = 2, main = "Averaged Periodogram from Bartlett's Method")
  list(combined_periodogram, maximumfreq, N)
}
