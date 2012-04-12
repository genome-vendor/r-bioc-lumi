\name{hist-methods}
\docType{methods}
\alias{hist-methods}
\alias{hist.ExpressionSet}
\alias{hist,ExpressionSet-method}

\title{Density plot of a ExpressionSet object}
\description{ Creating density plot of sample intensities in a ExpressionSet object. It is equivalent to \link{density-methods}. }

\usage{
\S4method{hist}{ExpressionSet}(x, ...)
}

\arguments{
	\item{x}{ a \code{\link[Biobase]{ExpressionSet}} object}
	\item{\dots}{ other parameters for \link{density-methods} function }
}

\seealso{
  \code{\link{LumiBatch-class}}, \link{density-methods}, \code{\link[graphics]{hist}}
}
\examples{
## load example data
data(example.lumi)

hist(example.lumi)
}
\keyword{methods}
\keyword{hplot}
