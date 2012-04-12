\name{lumiQ}
\alias{lumiQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Quality control evaluation of the LumiBatch object }
\description{
Quality control evaluation of the LumiBatch object and returns a summary of the data
}

\usage{
lumiQ(x.lumi, logMode = TRUE, detectionTh = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.lumi}{ a LumiBatch object }
  \item{logMode}{ transform as log2 or not (the function can check whether it is already log transformed.)}
  \item{detectionTh}{ the detection threshold used by \code{\link{detectionCall}} } 
}
\details{
  Quality control of a LumiBatch object includes estimating the mean and standard deviation of the chips, detectable probe ratio of each chip, sample (chip) relations, detecting outliers of samples (chips). The produced QC information is kept in the QC slot of LumiBatch class. The summary function will provide a summary of the QC information (See example).
}
\value{
a LumiBatch object with QC slot keeping the QC information
}

\author{ Pan Du }

\seealso{\code{\link{LumiBatch}}, \code{\link{plot,LumiBatch-method}} }
\examples{
## load example data
data(example.lumi)

## Do quality control estimation
lumi.Q <- lumiQ(example.lumi)

## A summary of the QC
summary(lumi.Q, 'QC')

## Plot the results
## plot the pairwise sample correlation
plot(lumi.Q, what='pair')

## see more examples in "plot,LumiBatch-method" help documents

}
\keyword{methods}
