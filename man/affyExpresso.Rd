\name{affyExpresso}
\alias{affyExpresso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Preprocess Affymetrix data by integrating VST with expresso method }
\description{
Preprocess Affymetrix data by integrating VST with expresso method
}
\usage{
affyExpresso(afbatch, bg.correct = TRUE, bgcorrect.method = NULL, bgcorrect.param = list(), variance.stabilize = TRUE, varianceStabilize.method = c("vst", "log2", "cubicRoot"), varianceStabilize.param = list(), normalize = TRUE, normalize.method = NULL, normalize.param = list(), pmcorrect.method = NULL, pmcorrect.param = list(), summary.method = NULL, summary.param = list(), summary.subset = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{afbatch}{ a vector of CEL file names or an AffyBatch object, see \code{\link[affy]{AffyBatch-class}}  }
  \item{bg.correct}{ a boolean to express whether background correction is wanted or not }
  \item{bgcorrect.method}{ the name of the background adjustment method }
  \item{bgcorrect.param}{ a list of parameters for bgcorrect.method (if needed/wanted) }
  \item{variance.stabilize}{ a boolean to express whether variance stabilization is wanted or not }
  \item{varianceStabilize.method}{ the name of the variance stabilizing transform, same as \code{\link{lumiT}} function }
  \item{varianceStabilize.param}{ a list of parameters for transformation method  }
  \item{normalize}{ normalization step wished or not }
  \item{normalize.method}{ the normalization method to use }
  \item{normalize.param}{ a list of parameters to be passed to the normalization method (if wanted) }
  \item{pmcorrect.method}{ the name of the PM adjustement method }
  \item{pmcorrect.param}{ a list of parameters for pmcorrect.method (if needed/wanted) }
  \item{summary.method}{ the method used for the computation of expression values }
  \item{summary.param}{ a list of parameters to be passed to the summary.method (if wanted) }
  \item{summary.subset}{ a list of 'affyids'. If NULL, a expression summary value is computed for everything on the chip }
  \item{verbose}{ logical value. If TRUE it writes out some messages }
}
\details{
 This function basically integrates the VST (variance stabilizing transformation) transformation into the expresso function in the affy package. The variance stabilization is based on the mean and variance relations of pixel intensities of each probe.
}
\value{
Return an object of class ExpressionSet. 
}
\references{  }
\author{ Pan Du }
\note{ 
The performance of this function is still under evaluation.
}
\seealso{ \code{\link[affy]{rma}} and \code{\link{vst}} }
\examples{


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
