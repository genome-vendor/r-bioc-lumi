\name{detectOutlier}
\alias{detectOutlier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Detect the outlier sample (or gene)}
\description{ 
Detect the outlier sample (or gene) based on distance to the cluster center
}

\usage{
detectOutlier(x, metric = "euclidean", standardize = TRUE, Th = 2, ifPlot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a LumiBatch object, ExpressionSet object or a matrix with each column corresponding to a sample or other profile }
  \item{metric}{ the distance matric}
  \item{standardize}{ standardize the profile or not }
  \item{Th}{ the threshold of outlier,  }
  \item{ifPlot}{ to plot the result (as a hierarchical tree) or not }
}
\details{
The current outlier detection is based on the distance from the sample to the center (average of all samples).
The assumption of the outlier detection is that there is only one single cluster and the distance from the sample
to the center is Gaussian distributed.

The outlier is detected when its distance to the center is larger than a certain threshold.
The threshold is calculated as Th * median distances to the center.

The profile relations can be visualized as a hierarchical tree.
}
\value{
Plot the results or	return the outlier (a logic vector) with the distance matrix and threshold as attributes.
}

\author{ Pan Du }

\seealso{ \code{\link{lumiQ}} }
\examples{
## load example data
data(example.lumi)

## detect the outlier (Further improvement needed.)
temp <- detectOutlier(example.lumi, ifPlot=TRUE)

}
\keyword{ methods }
