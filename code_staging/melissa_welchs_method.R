<<<<<<< Updated upstream
  welchs_method <- function(x, seglength, overlap = 0) {
    # Catches
    stopifnot(is.numeric(x), length(x) > 1)
    if (!is.ts(x)) {
      x <- as.ts(x)
      warning("x is not a time-series, conversion 'as.ts(x)' was applied.")
    }
    
    # splitWithOverlap function found on stackexchange - splits a vector into segments of length seg.length
    # with overlap of overlap
    splitWithOverlap <- function(vec, seg.length, overlap) {
      starts = seq(1, length(vec), by=seg.length-overlap)
      ends   = starts + seg.length - 1
      ends[ends > length(vec)] = length(vec)
      
      lapply(1:length(starts), function(i) vec[starts[i]:ends[i]])
    }
    
    # split the time series
    splitup <- splitWithOverlap(x, seglength, overlap)
    first_bad <- min(which(unlist(lapply(splitup, FUN = length)) < seglength))
    splitup <- splitup[1:first_bad]
    L <- length(splitup) # how many segments are there?
    splitup[[L]] <- x[(length(x)-seglength+1):length(x)]
    
    # create a Welch window function, whose definition was found on Wikipedia. 
    welch_window <- function(N) {
      n <- 1:N
      welch <- 1 - ((n - (N-1)/2)/((N-1)/2))^2
      return(welch)
    }
    
    # multiply each of the segments by the Welch Window function
    for (i in 1:L) {
      splitup[[i]] <- splitup[[i]]*welch_window(seglength)
    }
    
    # calculate the periodogram for each of the segments 
    M <- length(splitup[[1]]) # what is the length of each segment? basically redefine seglength
    zerop <- 2^(ceiling(log(M, 2)) + 1) - M
    nFFT <- zerop + M
    w <- (zerop+M)/2 + 1 # time to calculate the periodograms for each of the segments 
    freq <- (0:(nFFT/2))/nFFT
             
    # Zero-padding, re-phrasing
    splitup <- lapply(splitup, FUN = function(x) { c(x, rep(0, zerop)) })
             
    # P <- list()
    # meanP <- matrix(data = NA, nrow = w, ncol = L) 
    # meanvector <- vector(length = w)
    # for (k in 1:w) {
    #  for (i in 1:L) {
    #    dft <- abs(fft(splitup[[i]])/sqrt(N))^2
    #    P[[i]] <- (4/N)*dft[1:w]
    #    meanP[k,i] <- P[[i]][k]
    #    meanvector[k] <- mean(meanP[k,])
    #  } }
    S <- vector("list", length = L)   
    for(i in 1:L) {
      S[[i]] <- (1 / sqrt(M)) * (abs(fft(splitup[[i]]))^2)[1:w]
    }
    S <- Reduce("+", S) / L
             
    # finally, plot the average of all the periodograms
    combined_periodogram <- plot(freq, S, type = "l", ylab = "average of log(I(lambda))", 
                                 lwd = 2, main = "Averaged Periodogram from Welch's Method",
                                 log = 'y')
    
    # return a list
    list(freq = freq, spec = S, nseg = L)
=======
welchs_method <- function(x, seglength, overlap = 0) {
  # Catches
  stopifnot(is.numeric(x), length(x) > 1)
  if (!is.ts(x)) {
    x <- as.ts(x)
    warning("x is not a time-series, conversion 'as.ts(x)' was applied.")
  }
  
  # splitWithOverlap function found on stackexchange - splits a vector into segments of length seg.length
  # with overlap of overlap
  splitWithOverlap <- function(vec, seg.length, overlap) {
    starts = seq(1, length(vec), by=seg.length-overlap)
    ends   = starts + seg.length - 1
    ends[ends > length(vec)] = length(vec)
    
    lapply(1:length(starts), function(i) vec[starts[i]:ends[i]])
  }
  
  # split the time series
  splitup <- splitWithOverlap(x, seglength, overlap)
  first_bad <- min(which(unlist(lapply(splitup, FUN = length)) < seglength))
  splitup <- splitup[1:first_bad]
  L <- length(splitup) # how many segments are there?
  splitup[[L]] <- x[(length(x)-seglength+1):length(x)]
  
  # create a Welch window function, whose definition was found on Wikipedia. 
  welch_window <- function(N) {
    n <- 1:N
    welch <- 1 - ((n - (N-1)/2)/((N-1)/2))^2
    return(welch)
>>>>>>> Stashed changes
  }
  
  # multiply each of the segments by the Welch Window function
  for (i in 1:L) {
    splitup[[i]] <- splitup[[i]]*welch_window(seglength)
  }
  
  # calculate the periodogram for each of the segments 
  M <- length(splitup[[1]]) # what is the length of each segment? basically redefine seglength
  zerop <- 2^(ceiling(log(M, 2)) + 1) - M
  nFFT <- zerop + M
  w <- (zerop+M)/2 + 1 # time to calculate the periodograms for each of the segments 
  freq <- (0:(nFFT/2))/nFFT
  
  # Zero-padding, re-phrasing
  splitup <- lapply(splitup, FUN = function(x) { c(x, rep(0, zerop)) })
  
  # P <- list()
  # meanP <- matrix(data = NA, nrow = w, ncol = L) 
  # meanvector <- vector(length = w)
  # for (k in 1:w) {
  #  for (i in 1:L) {
  #    dft <- abs(fft(splitup[[i]])/sqrt(N))^2
  #    P[[i]] <- (4/N)*dft[1:w]
  #    meanP[k,i] <- P[[i]][k]
  #    meanvector[k] <- mean(meanP[k,])
  #  } }
  S <- vector("list", length = L)   
  for(i in 1:L) {
    S[[i]] <- (1 / sqrt(M)) * (abs(fft(splitup[[i]]))^2)[1:w]
  }
  S <- Reduce("+", S) / L
  
  # finally, plot the average of all the periodograms
  combined_periodogram <- plot(freq, S, type = "l", ylab = "average of log(I(lambda))", 
                               lwd = 2, main = "Averaged Periodogram from Welch's Method",
                               log = 'y')
  
  # return a list
  list(freq = freq, spec = S, nseg = L)
}