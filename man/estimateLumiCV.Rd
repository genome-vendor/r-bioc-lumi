\name{estimateLumiCV}
\alias{estimateLumiCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimate the coefficient of variance matrix of LumiBatch object }
\description{ 
Estimate the coefficient of variance matrix of LumiBatch object for each measurement or probe.
}

\usage{
estimateLumiCV(x.lumi, type = c("measurement", "probe"), ifPlot = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.lumi}{ a LumiBatch object }
  \item{type}{ estimate the coefficient of variance of each measurement or each probe }
  \item{ifPlot}{ determince whether to plot the density plot or not }
  \item{\dots}{ optional arguments to \code{\link[graphics]{plot}}.}
}
\details{
  By default, the coefficient of variance is the ratio of the mean and variance of the bead expression values.
Basically, it is the ration of \code{exprs}  and \code{se.exprs} element  of \code{LumiBatch} object.
If the type is "probe", it is the ratio of the mean and variance of probe expression profile.
}
\value{
A matrix of coefficient of variance
}

\author{Pan Du}

\seealso{ \code{\link{lumiQ}}}
\examples{
## load example data
data(example.lumi)

## estimate the coefficient of variance and plot the density plot of it
cv <- estimateLumiCV(example.lumi, ifPlot = TRUE)

}
\keyword{ methods }
