\name{lumiN}
\alias{lumiN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Between chip normalization of a LumiBatch object }
\description{
  A main function of between chip normalization of a LumiBatch object. Currently, four methods ("rsn", "ssn", "quantile", "loess", "vsn") are supported.
}
\usage{
lumiN(x.lumi, method = c("quantile", "rsn", "ssn", "loess", "vsn"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{x.lumi}{ an ExpressionSet inherited object or a data matrix with columns as samples and rows as genes}
  	\item{method}{ four different between chips normalization methods ("quantile", "rsn", "ssn", "loess", "vsn") are supported }
  	\item{\dots}{ other parameters used by corresponding method }
}
\details{
lumiN is an interface for different normalization methods. Currently it supports "RSN" (See \code{\link{rsn}}), "SSN" (See \code{\link{ssn}}), "loess" (See \code{\link[affy]{normalize.loess}}), "quantile" (See \code{\link[affy]{normalize.quantiles}}) and "VSN" (See \code{\link[vsn]{vsn}}). See details in individual functions. Note: the "VSN" normalization should be directly applied to the raw data instead of the lumiT processed data.
}

\value{
Return an object with expression values normalized. The class of the return object is the same as the input object x.lumi.
If it is a LumiBatch object, it also includes the VST transform function and its parameters as attributes: "transformFun", "parameter". See \code{\link{inverseVST}} for details.
}
\references{  }
\author{Pan Du, Simon Lin}
\note{ 

}
\seealso{  \code{\link{rsn}}, \code{\link{ssn}} }
\examples{
## load example data
data(example.lumi)

## Do lumi transform
lumi.T <- lumiT(example.lumi)

## Do lumi between chip normaliazation
lumi.N <- lumiN(lumi.T, method='rsn', ifPlot=TRUE)

}
\keyword{ methods }
