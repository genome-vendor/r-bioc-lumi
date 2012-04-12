\name{affyVstRma}
\alias{affyVstRma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Preprocess Affymetrix data by integrating VST with RMA method }
\description{
Preprocess Affymetrix data by integrating VST with RMA method
}
\usage{
affyVstRma(afbatch, bgcorrect.method = "none", bgcorrect.param = list(), VST.param = list(), verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{afbatch}{ a vector of CEL file names or an AffyBatch object, see \code{\link[affy]{AffyBatch-class}} }
  \item{bgcorrect.method}{ the name of the background adjustment method }
  \item{bgcorrect.param}{ a list of parameters for bgcorrect.method (if needed/wanted) }
  \item{VST.param}{ a list of parameters for \code{\link{vst}} method }
  \item{verbose}{ logical value. If TRUE it writes out some messages. }
  \item{\dots}{ other parameters used by \code{\link[affy]{rma}} function }
}
\details{
This function basically integrates the VST (variance stabilizing transformation) transformation into the rma function in the affy package. The variance stabilization is based on the mean and variance relations of pixel intensities of each probe.
}
\value{
 Return an object of class ExpressionSet. 
}
\references{  }
\author{ Pan Du }
\note{ 
The performance of this function is still under evaluation.
}
\seealso{ \code{\link[affy]{expresso}} and \code{\link{vst}} }
\examples{


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
