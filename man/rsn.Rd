\name{rsn}
\alias{rsn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Robust Spline Normalization between chips  }
\description{
  Robust spline normalization (monotonic curves) between chips
}
\usage{
rsn(x.lumi, targetArray = NULL, excludeFold = 2, span = 0.03, ifPlot = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.lumi}{ an ExpressionSet inherited object or a data matrix with columns as samples and rows as genes }
  \item{targetArray}{ A target chip is the model for other chips to normalize. It can be a column index, a vector or a LumiBatch object with one sample. }
  \item{excludeFold}{ exclude the genes with fold change larger than "excludeFold" during fitting the curve in normalization }
  \item{span}{ the span parameter used by \code{\link{monoSmu}} }
  \item{ifPlot}{ determine whether to plot intermediate results }
  \item{\dots} { other parameters used by \code{\link{monoSmu}} }
}
\details{
The robust spline normalization (RSN) algorithm combines the features of quantile and loess normalization. It is designed to normalize the variance-stabilized data. The function will check whether the data is variance stabilized (vst or log2 transform), if not, it will automatically run lumiT before run rsn. For details of the algorithm, please see the reference.

The targetArray can be a column index, a vector or a LumiBatch object with one sample, which corresponds to an external sample to be normalized with. This is very useful for handling large data set or normalizing the data set with a common reference (targetArray).
}
\value{
Return an object with expression values normalized. The class of the return object is the same as the input object x.lumi.
If it is a LumiBatch object, it also includes the VST transform function and its parameters as attributes: "transformFun", "parameter". See \code{\link{inverseVST}} for details.
}
\references{ 

}
\author{Pan Du, Simon Lin}
\note{ 

}
\seealso{ \code{\link{lumiN}}, \code{\link{monoSmu}} }
\examples{

}
\keyword{ methods }
