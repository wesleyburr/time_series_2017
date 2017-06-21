\name{acf_arma}
\alias{acf_arma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Autocorrelation function estimation for casual ARMA(p,q) process: acf_arma()
}
\description{
The function acf_arma() computes estimates of the autocorrelation/autocoveriance up to given lag.
It retruns a data frame contain lag (lags), autocovariance (gamma) and autocorrelation (rho). 
}
\usage{
acf_arma <- function( phi, theta, lag = 20, dig = 5 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phi}{A vector, containing coefficient of AR part as array. it require to have 1 as leading coefficient}
  \item{theta}{A vector, containing coefficient of MA part as array. it require to have 1 as leading coefficient}	
  \item{lag}{the desired maximun lag , The default value is 20.}
  \item{dig}{the desired rounding digit of return value, the default value is 5}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
If an integer (within 8 decimal places) is not entered for the value of dig,
then the floor of dig will be used instead.
If an positive integer (within 8 decimal places) is not entered for the value of lag,
then the floor of absolute value of lag will be used instead.

\value{
\item{lags} {Corresponding lag}
\item{gamma} {Autocovariance estimation}
\item{rho} {Autocorrelation estimation}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Brockwell, P. J., & Davis, R. A. (2013). Time series: theory and methods. New York, NY: Springer.
%% ~put references to the literature/web site here ~
}
\author{
Ankai Liu
}
Melissa Van Bussel, Carlone Scott, Skye Griffith, Steph Johnson, Mark Weygang, Ankai Liu, Sam Leigh, Wesley Burr. 
\note{
%%  ~~further notes~~
}
\Warning{
%% ~Make other sections like Warning with \section{Warning }{....} ~
xd
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~

}
\examples{
# Define time series data in a vector, x
phi <- c(1, -1, 1 / 4)
theta <- c(1, 1)

# Use acf_arma function 
acf_arma(phi, theta)
acf_arma(phi, theta, 20, 8)

# Another example
phi <- c(1, rnorm(100))
theta <- c(1,rnorm(50))
 
acf_arma(phi, theta)
acf_arma(phi, theta, 10, 8)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~autocorrelation }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~autocovariance }
\keyword{ ~acf }% __ONLY ONE__ keyword per line
\keyword{ ~ARMA }