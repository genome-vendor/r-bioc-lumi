\name{monoSmu}
\alias{monoSmu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Monotonic smooth method }
\description{
  Fit the monotonic-constraint spline curve
}
\usage{
monoSmu(x, y, newX = NULL, nSupport = min(200, length(x)), nKnots = 6, rotate = FALSE, ifPlot = FALSE, xlab = 'x', ylab = 'y', ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a vector represents x values }
  \item{y}{ a vector represents y values }
  \item{newX}{ the new values to be transformed. If not provided, "x" will be used. }
  \item{nSupport}{ downsampled data points }
  \item{nKnots}{ parameter used by \code{\link{monoSpline}} }
  \item{rotate}{ determine whether to rotate the axis with 45 degrees in clockwise, i.e., fit the curve in the MA-plot. }
  \item{ifPlot}{ determine whether to plot intermediate results }
  \item{xlab}{ the xlab of the plot }
  \item{ylab}{ the ylab of the plot }
  \item{\dots}{ parameters used by \code{\link{supsmu}} and \code{\link{plot}} }
}
\details{
function called by lumiN.rsn. The function first fits a monotonic spline between vector x and y, then transforms the vector newX based on the fitted spline. (After transformation the fitted spline is supposed to be a diagonal line, i.e., x=y)
}
\value{
Return the transformed "newX" based on the smoothed curve
}
\references{ 
Lin, S.M., Du, P., Kibbe, W.A.,  {\it Model-based Variance-stabilizing Transformation for Illumina Microarray Data}, submitted
}
\author{ Simon Lin, Pan Du }

\seealso{ \code{\link{monoSpline}}  }
\examples{

x <- sort(runif(200) * 100)
y <- x + rnorm(200) * 10
temp <- monoSmu(x, y, ifPlot=TRUE)

}
\keyword{ methods }
