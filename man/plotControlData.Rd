\name{plotControlData}
\alias{plotControlData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot the mean expression (with standard deviation bar) of different type of control probes }
\description{
  Plot the mean expression (with standard deviation bar) of different type of control probes. Multiple control types can be plotted in a single plot. The available control types can be get by running getControlType(controlData).
}
\usage{
plotControlData(controlData, type = NULL, slideIndex = NULL, logMode = FALSE, new = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{controlData}{ a LumiBatch object including control data or a control data data.frame}
  \item{type}{ the control probe type (case insensitive), which can be get by running getControlType(controlData) }
  \item{slideIndex}{ the slide index or ID corresponding to each sample }
  \item{logMode}{ whether show the data in log2 scale }
  \item{new}{ whether refresh the new plot or add it on the old one  }
  \item{\dots}{ other parameters used by default plot function }
}
\details{
  When multiple control types are selected, they will be plotted in a two-column plot.
}
\value{
	plot the picture and return TRUE if everything is OK
}
\references{  }
\author{ Pan Du }
\seealso{ \code{\link{addControlData2lumi}} }
\examples{
	controlFile <- system.file('doc', 'Control_Probe_Profile.txt', package='lumi')
	controlData <- getControlData(controlFile)
	getControlType(controlData)
	plotControlData(controlData, type='NEGATIVE')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{hplot}
