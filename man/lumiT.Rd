\name{lumiT}
\alias{lumiT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Transfer the Illumina data to stabilize the variance }
\description{
  Transfer the Illumina data to stabilize the variance.
}
\usage{
lumiT(x.lumi, method = c("vst", 'log2', 'cubicRoot'), ifPlot = FALSE, stdCorrection = TRUE, simpleOutput = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.lumi}{ LumiBatch object }
  \item{method}{ four methods are supported: "vst", "log2", "cubicRoot" }
  \item{ifPlot}{ determine whether to plot the intermediate results }
  \item{stdCorrection}{ determine transfer the standard error of the mean as the standard deviation, used for 'vst' method. }
  \item{simpleOutput}{ determine whether to simplify the output LumiBatch object, which will set the se.exprs, detection and beadNum slots as NULL. }
  \item{\dots}{ other parameters used by \code{\link{vst}} }
}
\details{
lumiT is an interface of difference variance stabilizing transformation. See \code{\link{vst}} for details of VST (Variance Stabilizing Transform) of Illumina data.

The adding of the parameter "stdCorrection" is for the value correction of the STDEV (or STDERR) columns when 'vst' method is selected. The STDEV (or STDERR) columns of the BeadStudio output file is the standard error of the mean of the bead intensities corresponding to the same probe. (Thanks Gordon Smyth kindly provided this information.). As the variance stabilization (see \code{\link{vst}} function) requires the information of the standard deviation instead of the standard error of the mean, the value correction is required. The corrected value will be x * sqrt(N), where x is the old value (standard error of the mean), N is the number of beads corresponding to the probe.
}

\value{
Return a LumiBatch object with transformed expression values. It also includes the VST transform function and its parameters as attributes: "transformFun", "parameter". See \code{\link{inverseVST}} for details.
	
}
\references{
 Lin, S.M., Du, P., Kibbe, W.A.,  {\it Model-based Variance-stabilizing Transformation for Illumina Microarray Data}, submitted
 }
\author{ Pan Du, Simon Lin }

\seealso{ \code{\link{vst}} }
\examples{

## load example data
data(example.lumi)

## Do default VST variance stabilizing transform
lumi.T <- lumiT(example.lumi, ifPlot=TRUE)


}
\keyword{ methods }
