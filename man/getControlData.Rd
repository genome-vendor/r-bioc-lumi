\name{getControlData}
\alias{getControlData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get control probe information }
\description{
  Get control probe information from Bead Studio output or a LumiBatch object.
}
\usage{
getControlData(x, type = c('data.frame', 'LumiBatch'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ the control data can be a LumiBatch object or the Control Probe Profile file outputted by BeadStudio }
  \item{type}{ determine the return data type }
}
\details{
  
}
\value{
   By default, it returns a data.frame with first two columns as "controlType" and "ProbeID". The rest columns are the expression amplitudes for individual samples. When type is 'LumiBatch', it returns a LumiBatch object, which basically is the return of lumiR without combining duplicated TargetIDs. As the return is a LumiBatch object, it includes more information, like probe number, detection p-value and standard error of the measurement.
}
\references{  }
\author{ Pan Du }
\seealso{ \code{\link{addControlData2lumi}} }
\examples{
	controlFile <- system.file('doc', 'Control_Probe_Profile.txt', package='lumi')
	## return a data.frame
	controlData <- getControlData(controlFile)
	class(controlData)
	names(controlData)
	
	## return a LumiBatch object
	controlData <- getControlData(controlFile, type='LumiBatch')
	summary(controlData)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
