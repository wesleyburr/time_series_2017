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
    splitup <- splitup[1:(length(splitup)-1)] # for some reason the splitWithOverlap function creates one too many segments
    L <- length(splitup) # how many segments are there?
    
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
    freq <- (0:(M/2))/M
    zerop <- 2^(ceiling(log(M, 2)) + 1) - M
    w <- (zerop+M)/2 + 1 # time to calculate the periodograms for each of the segments 
             
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
      S[[i]] <- (4/M) * (abs(fft(splitup[[i]])/sqrt(M))^2)[1:w]
    }
    S <- Reduce("+", S) / L
             
    # finally, plot the average of all the periodograms
    combined_periodogram <- plot(freq, S, type = "h", ylab = "average of log(I(lambda))", 
                                 lwd = 2, main = "Averaged Periodogram from Welch's Method"
                                 log = 'y')
    
    # return a list
    list(combined_periodogram, L)
  }
  
