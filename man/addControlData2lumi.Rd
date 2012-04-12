\name{addControlData2lumi}
\alias{addControlData2lumi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Add the control probe data into the controlData slot of LumiBatch object }
\description{
   Add the control probe profile data, outputted by BeadStudio, into the controlData slot of LumiBatch object. 
}
\usage{
addControlData2lumi(controlData, x.lumi)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{controlData}{ the control data can be a data.frame or the control probe filename outputted by BeadStudio }
  \item{x.lumi}{ a LumiBatch object, to which controlData will be added. }
}
\details{
 The controlData slot in LumiBatch object is a data.frame with first two columns as "controlType" and "ProbeID". The rest columns are the expression amplitudes for individual samples.
}
\value{
Return the LumiBatch object with controlData slot filled.
}
\references{  }
\author{ Pan Du }

\seealso{  \code{\link{getControlData}}, \code{\link{plotControlData}} }
\examples{
	## Not runnable
	# controlFile <- 'Control_Probe_Profile.txt'
	# x.lumi <- addControlData2lumi(controlFile, x.lumi)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
