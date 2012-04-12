\name{plotStringencyGene}
\alias{plotStringencyGene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ plot the Stringency related control probe profiles }
\description{
  Plot the Stringency related control probe (Low-Stringency, Medium-Stringency and High-Stringency) profiles. Using getControlType function to view available stringency types.
}
\usage{
plotStringencyGene(controlData, lib = NULL, slideIndex = NULL, addLegend = TRUE, logMode = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{controlData}{ a LumiBatch object including control data or a control data data.frame }
  \item{lib}{ the annotation library (for retrieving the gene name) }
  \item{slideIndex}{ the slide index or ID corresponding to each sample }
  \item{addLegend}{ whether add legend or not }
  \item{logMode}{ whether show the data in log2 scale  }
  \item{\dots}{ other parameters used by default matplot function }
}
\details{
  
}
\value{
	plot the picture and return TRUE if everything is OK
}
\references{  }
\author{ Pan Du }
\seealso{ \code{\link{addControlData2lumi}}, \code{\link{plotControlData}} }

\examples{
	controlFile <- system.file('doc', 'Control_Probe_Profile.txt', package='lumi')
	controlData <- getControlData(controlFile)
	plotStringencyGene(controlData)
}
\keyword{methods}
\keyword{hplot}
