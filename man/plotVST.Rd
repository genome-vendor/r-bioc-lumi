\name{plotVST}
\alias{plotVST}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ plot the VST (Variance Stabilizing Transform) function }
\description{
plot the VST (Variance Stabilizing Transform) function of VST transformed LumiBatch object or parameters of VST function.
}
\usage{
plotVST(x, transFun = NULL, plotRange = NULL, addLegend = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a LumiBatch object after lumiT transform, or a matrix or data.frame with VST parameter }
  \item{transFun}{ a character vector of transformation function (asinh or log2) }
  \item{plotRange}{ the plot range of untransformed data }
  \item{addLegend}{ add legend or not }
  \item{\dots}{ other parameter used by \code{\link[graphics]{plot}} function  }
}
\details{
  

}
\value{
invisibly return the untransformed and transformed values.
}
\author{ Pan Du }
\seealso{  \code{\link{vst}} }
\examples{

## load example data
data(example.lumi)

## Do default VST variance stabilizing transform
lumi.T <- lumiT(example.lumi, ifPlot=TRUE)

## plot the transform function
plotVST(lumi.T)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
