\name{lumiExpresso}
\alias{lumiExpresso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{From raw Illumina probe intensities to expression values}
\description{
  Goes from raw Illumina probe intensities to expression values
}
\usage{
lumiExpresso(lumiBatch, bg.correct = FALSE, bgcorrect.param = list(), variance.stabilize = TRUE, 
	varianceStabilize.param = list(), normalize = TRUE, normalize.param = list(), QC.evaluation = TRUE, 
	QC.param = list(), verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lumiBatch}{ a LumiBatch object, which can be the return of \code{\link{lumiR}} }
  \item{bg.correct}{ a boolean to decide whether to do background correction or not }
  \item{bgcorrect.param}{ a list of parameters of \code{\link{lumiB}} }
  \item{variance.stabilize}{ a boolean to decide whether to do variance stabilization or not }
  \item{varianceStabilize.param}{ a list of parameters of \code{\link{lumiT}} }
  \item{normalize}{ a boolean to decide whether to do normalization or not }
  \item{normalize.param}{ a list of parameters of \code{\link{lumiN}} }
  \item{QC.evaluation}{ a boolean to decide whether to do quality control estimation before and after preprocessing }
  \item{QC.param}{ a list of parameters of \code{\link{lumiQ}} }
  \item{verbose}{ a boolean to decide whether to print out some messages }
}
\details{
  The function is to encapsulate the major functions of Illumina preprocessing. It is organized in a similar way as the \code{\link[affy]{expresso}} function in affy package.
}
\value{
  return a processed LumiBatch object. The operation history can be track in the history slot of the object.
}

\author{ Pan Du }

\seealso{ \code{\link{lumiB}}, \code{\link{lumiT}}, \code{\link{lumiN}} }
\examples{
## load example data
data(example.lumi)

## Do all the default preprocessing in one step
lumi.N <- lumiExpresso(example.lumi)

## Do customized preprocessing. No variance stabilizing or log transform, use Quantile normalization.
lumi.N <- lumiExpresso(example.lumi, variance.stabilize=FALSE, normalize.param = list(method='quantile'))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
