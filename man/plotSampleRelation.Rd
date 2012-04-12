\name{plotSampleRelation}
\alias{plotSampleRelation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ visualize the sample relations }
\description{
  plot the sample relations based on MDS or hierarchical clustering
}
\usage{
plotSampleRelation(x, selProbe = NULL, cv.Th = 0.1, standardize = TRUE, method = c("cluster", "mds"), dimension = c(1, 2), color = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a LumiBatch object, ExpressionSet object or a matrix with each column corresponding to a sample}
  \item{selProbe}{ the selected probes used to determine the sample relations. If not provide, all the probes will be used. }
  \item{cv.Th}{ the threshold of the coefficient of variance of probes used to select probes to estimate sample relations }
  \item{standardize}{ standardize the expression profiles or not }
  \item{method}{ "MDS" or "hierarchical clustering" }
  \item{dimension}{ the principle components to visualize the MDS plot }
  \item{color}{ the color for each sample during plot. Only support the "mds" method }
}
\details{
  Estimate the sample relations based on selected probes (based on large coefficient of variance (mean / standard variance)). Two methods can be used: MDS (Multi-Dimensional Scaling) or hierarchical clustering methods.
}
\value{
Plot the results or return the distance matrix.
}

\author{ Pan Du }

\seealso{  \code{\link{lumiQ}}, \code{\link{LumiBatch}}, , \code{\link{plot.LumiBatch}} }
\examples{

## load example data
data(example.lumi)

## plot the sample relations with MDS
## the color of sample is automatically set based on the sample type
plotSampleRelation(example.lumi, col=c('100US', '95US:5P', '100US', '95US:5P'))

## plot the sample relations with hierarchical clustering
plotSampleRelation(example.lumi, method='cluster')

}
\keyword{ hplot }
