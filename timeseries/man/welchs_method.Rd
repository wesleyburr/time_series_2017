\name{welchs_method}
\alias{welchs_method}
\title{Welch's Method Calculator}
\usage{
welchs_method(x, seglength, overlap)
}
\arguments{
\item{x}{A time series (numeric vector)}
\item{seglength}{The length of each segment that the time series should be split up into.}
\item{overlap}{The number of points that each segment should overlap by.}
}
\value{
The averaged periodogram computed using Welch's Method, as well as the maximum value of the periodogram, and the number of segments.
}
\description{
Takes a time series vector, splits it up into segments, and computes an averaged periodogram via Welch's Method.
}
\note{
Uses a Welch Window function.
}
\examples{
data(sunspots)
dat <- sunspots
segmentlength <- 282
d <- 141 # 50 percent overlap
welchs_method(x = dat, seglength = segmentlength, overlap = d)
}
\references{
Brockwell, P. J., & Davis, R. A. (2013). Time series: Theory and Methods. New York, NY: Springer.
}
\author{
Melissa Van Bussel <melissavanbussel@trentu.ca>
}
\keyword{Welch}
\keyword{Method}
\keyword{Spectral Analysis}
\keyword{Periodogram}
\keyword{Window Function}