fftperiodogram <- function(x, xfreq = frequency(x)) {
  # Catches and warnings 
  stopifnot(is.numeric(x), length(x) > 0, is.numeric(xfreq))
  
  name <- deparse(substitute(x))
  n <- length(x)
  xfreq <- frequency(x)
  Nspec <- floor(n/2)
  freq <- seq.int(from = xfreq/n, by = xfreq/n, length.out = Nspec)
 
  # Zero-pad the series 
  zeropad <- 2^(ceiling(log(N, 2)) + 1) - N
  x <- c(x, rep(0, zeropad))
  
  # Reassign length of series 
  N <- length(x)
  w <- N/2 + 1
  
  # Calculate periodogram and plot it
  dft <- (1/sqrt(N))*(abs(fft(x))^2)[1:w]
  dft <- log(dft)
  periodogramplot <- plot(freq, dft[1:Nspec], type = "l", lwd = 1, main = paste("Periodogram of", name, "dataset", sep = " "), 
                          xlab = "Frequency",
                          ylab = "log of periodogram",
                          cex.main = 0.8)
  
  # Find frequency at which maximum value occurs
  freq_which_max <- freq[max(dft[2:length(dft)])]
  periodogram_return <- list(plot = periodogramplot, freq_which_max = freq_which_max)
  return(periodogram_return)
}


