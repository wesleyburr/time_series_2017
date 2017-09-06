\name{directspec}
\alias{directspec}
\title{Direct Spectrum Periodogram Calculator}
\usage{
directspec(x, xfreq = frequency(x), window = "welch", lwd = 1)
}
\arguments{
\item{x}{A time series (numeric vector)}
\item{xfreq}{The frequency of the time series; defaults to frequency(x). Changes the scale on the x-axis.}
\item{window}{The window function to be used for the calculation. Options are "welch", "hann", "hamming", "nutall". Defaults to "welch".}
\item{lwd}{The line width of the plot (must be numeric). Defaults to 1.}
}
\value{
A list which outputs the following:
\item{freq}{The frequencies used in the periodogram}
\item{plot}{The plotted periodogram}
\item{spec}{The values on the y-axis of the plotted periodogram for each of the frequencies in freq. These values are from
multiplying the time series by the specified window function, zero-padding, computing the square of the Discrete Fourier Transform, and then taking the log.}
\item{freq_which_max}{The frequency at which the maximum value in the computed periodogram occurred.}
}
\description{
Takes a time series vector and computes a direct spectrum estimate (periodogram): multiplies the series by a specified window function, zero-pads the result, finds the square of the Discrete Fourier Transform of this, and then takes the log. This is then plotted.
}
\note{
The xfreq argument is especially useful when the data you have is not of ts format. 
}
\examples{
data(sunspots)
dat <- sunspots
windowtype <- "hann" 
linewidth <- 2 
directspec(x = dat, window = windowtype, lwd = linewidth)
}
\references{
Brockwell, P. J., & Davis, R. A. (2013). Time series: Theory and Methods. New York, NY: Springer.
}
\author{
Melissa Van Bussel <melissavanbussel@trentu.ca>
}
\keyword{Periodogram}
\keyword{Direct Spectrum Estimation}
\keyword{Spectral Analysis}
\keyword{Bartlett}
\keyword{Window Function}
\keyword{Welch}