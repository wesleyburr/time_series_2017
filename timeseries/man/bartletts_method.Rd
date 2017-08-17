\name{bartletts_method}
\alias{bartletts_method}
\title{Bartlett's Method Calculator}
\usage{
bartletts_method(x, segments)
}
\arguments{
\item{x}{A time series (numeric vector)}
\item{segments}{The number of segments that the time series, x, should be split up into}
}
\value{
The averaged periodogram computed using Bartlett's Method, as well as the maximum value of the periodogram, and the length of each segment.
}
\description{
Takes a time series vector, splits it up into segments, and computes an averaged periodogram via Bartlett's Method.
}
\note{
If the length of the time series is not a multiple of the number of segments, zeroes will be added at the end of the time series.
}
\examples{
data(sunspot.month)
dat <- sunspot.month
num_segs <- 10
bartletts_method(x = dat, segments = num_segs)
}
\references{
Brockwell, P. J., & Davis, R. A. (2013). Time series: Theory and Methods. New York, NY: Springer.
}
\author{
Melissa Van Bussel <melissavanbussel@trentu.ca>
}
\keyword{Bartlett}
\keyword{Method}
\keyword{Spectral Analysis}
\keyword{Periodogram}
