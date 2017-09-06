directspec <- function(x, xfreq = frequency(x), window = "welch", lwd = 1) {
  # Catches
  stopifnot(is.numeric(x), length(x) > 1)
  name <- deparse(substitute(x))
  if (!is.ts(x)) {
    x <- as.ts(x)
    name <- "x"
    warning("x is not a time-series, conversion 'as.ts(x)' was applied.")
  }
  if (window != "welch" & window != "hann" & window != "hamming" & window != "nutall") {
    window <- "welch"
    warning("the window you entered was not valid, changing to a welch window instead")
  }
  if (!is.numeric(lwd)) {
    lwd <- 1
    warning("line width must be numeric, defaulted to 1 instead")
  }
  
  # length of time series and frequency 
  N <- length(x)
  Nspec <- floor(N/2)
  freq <- seq.int(from = xfreq/N, by = xfreq/N, length.out = Nspec)
  
  # Define a few different window functions 
  # Welch's window function
  welch_window <- function(N) {
    n <- 1:N
    welch <- 1 - ((n - (N-1)/2)/((N-1)/2))^2
    return(welch)
  }
  # Hann window function
  hann_window <- function(N) {
    n <- 1:n
    hann <- sin((pi*n) / (N-1))^2
    return(hann)
  }
  # Hamming window function
  hamming_window <- function(N) {
    n <- 1:N
    alpha <- 0.53836
    beta <- 0.46164
    hamming <- alpha - beta * cos((2*pi*n) / (N-1))
    return(hamming)
  }
  # Nutall window function 
  nutall_window <- function(N) {
    n <- 1:N
    a0 <- 0.355768
    a1 <- 0.487396
    a2 <- 0.144232
    a3 <- 0.012604
    nutall <- a0 - a1 * cos((2*pi*n) / (N-1)) + a2 * cos((4*pi*n) / (N-1)) - a3 * cos((6*pi*n) / (N-1))
    return(nutall)
  }
  # Multiply series by Welch's window function
  if (window == "welch") {
    windowtype <- welch_window(N)
  } else {
    if (window == "hamming") {
      windowtype <- hamming_window(N)
    } else {
      if (window == "hann") {
        windowtype <- hann_window(N)
      } else {
        if (window == "nutall") {
          windowtype <- nutall_window(N)
        }
      }
    }
  }
  x <- x*windowtype
  
  # Zero-pad the series 
  zeropad <- 2^(ceiling(log(N, 2)) + 1) - N
  x <- c(x, rep(0, zeropad))
  
  # Reassign length of series 
  N <- length(x)
  w <- N/2 + 1
  
  # Calculate periodogram and plot it
  dft <- (1/sqrt(N))*(abs(fft(x))^2)[1:w]
  dft <- log(dft)
  periodogramplot <- plot(freq, dft[1:Nspec], type = "l", lwd = lwd, main = paste("Periodogram of", name, "dataset using", window, "Window", sep = " "), 
       xlab = "Frequency",
       ylab = "log of periodogram",
       cex.main = 0.8)
  
  # Find frequency at which maximum value occurs
  freq_which_max <- freq[max(dft[2:length(dft)])]

  # Return a list
  list(freq = freq, plot = periodogramplot, spec = dft, freq_which_max = freq_which_max)
}