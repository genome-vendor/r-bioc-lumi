\name{monoSpline}
\alias{monoSpline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fitting a curve with monotonic spline }
\description{
Fitting a curve with monotonic spline
}
\usage{
monoSpline(x, y, newX=NULL, nKnots = 6, ifPlot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{x}{ a vector represents x values }
	\item{y}{ a vector represents y values }
	\item{newX}{ the new values to be transformed. If not provided, "x" will be used. }
	\item{nKnots}{ parameter used by function smoothCon in package mgcv }
  	\item{ifPlot}{ determine whether to plot intermediate results }
}
\details{
Function internally called by monoSmu
}
\value{
return the transformed "newX" based on the smoothed curve
}
\references{ }
\author{ Simon Lin, Pan Du }

\seealso{ \code{\link{monoSmu}}  }
\examples{

x <- sort(runif(200) * 100)
y <- x + rnorm(200) * 10
temp <- monoSpline(x, y, newX=1:100, ifPlot=TRUE)


}
\keyword{ methods }
