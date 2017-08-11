\name{fftperiodogram}
\alias{fftperiodogram}
\title{Periodogram Calculator}
\usage{
fftperiodogram(x, seconds)
}
\arguments{
\item{x}{A numeric vector of time series points, beginning at x_1 and ending at x_n}
\item{seconds}{The number of seconds between each observation in the time series; e.g., monthly data has a seconds value of 86400*30}
}
\value{
The periodogram for x, as well as the period which corresponds to the frequency with the highest periodogram value. 
}
\description{
Takes a vector of time series data, computes the periodogram via the discrete fourier transform (using the Fast Fourier Transform (fft) algorithm), and prints out the period which corresponds to the frequency with the highest periodogram value.
}
\note{
If the length of the time series is not an even number, the function will remove the first observation in the time series. There must be an even number of observations in order for the Fast Fourier Transform algorithm to work as we want it to. 
}
\examples{
# define the time series you'll be using
x_vector <- sunspot.year
num_seconds <- 86400*365

# Use function
fftperiodogram(x = x_vector, seconds = num_seconds)

}
\references{
Brockwell, P. J., & Davis, R. A. (2013). Time series: Theory and Methods. New York, NY: Springer.
}
\author{
Melissa Van Bussel <melissavanbussel@trentu.ca>
}
\keyword{Periodogram}
\keyword{Spectral Density}
\keyword{Spectral Analysis}
\keyword{FFT}
