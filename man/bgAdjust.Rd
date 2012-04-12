\name{bgAdjust}
\alias{bgAdjust}
\title{ Background adjustment for Illumina data }
\description{
 The method adjusts the data by subtracting an offset, which is estimated based on the quantile of the control probes
}
\usage{
bgAdjust(lumiBatch, probs = 0.5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lumiBatch}{ A LumiBatch object with controlData slot include control probe information }
  \item{probs}{ The quantile used to estimate the background }
  \item{\dots}{ other parameters used by \code{\link{quantile}} method }
}
\details{
The method adjusts the data by subtracting an offset, which is estimated based on the quantile of the control probes. The control probe information is kept in the controlData slot of the LumiBatch object. If no control data information, the method will do nothing.
}
\value{
  It returns a LumiBatch object with background adjusted.
}

\author{ Pan Du }

\seealso{ \code{\link{lumiB}} }
\examples{
	data(example.lumi)
	## Here will assume the minimum of the control probe as the background,
	##    because there is no negative control (blank beads) information for the Barnes data.
	example.lumi.b <- bgAdjust(example.lumi, probs=0)
}
\keyword{methods}

