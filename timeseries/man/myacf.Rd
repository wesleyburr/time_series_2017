\name{myacf}
\alias{myacf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Autocorrelation function estimation: myacf()
}
\description{
The function myacf computes estimates of the autocorrelation function. 
}
\usage{
myacf(x, maxLag = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector, containing a time series with at least two points.}
  \item{maxLag}{The highest lag for which you would like to calculate the acf. This must be an integer less than or equal to one minus the length of the time series vector, x. The default value is 20.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
If an integer (within 8 decimal places) is not entered for the value of maxLag, then the floor of maxLag will be used instead.
The myacf function differs from R's built in acf function because it does not plot the ACF and it does not compute the ACV. 
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
Melissa Van Bussel, Carlone Scott, Skye Griffith, Steph Johnson, Mark Weygang, Ankai Liu, Sam Leigh, Wesley Burr. 
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~

}
\examples{
# Define time series data in a vector, x
x <- seq(1:30)

# Use myacf function 
acf(x, maxLag = 12)
acf(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~correlation }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~acf }% __ONLY ONE__ keyword per line
