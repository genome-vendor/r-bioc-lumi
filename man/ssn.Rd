\name{ssn}
\alias{ssn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simple Scaling Normalization }
\description{
  This function basically adjusts the samples to the same background level and then optionally scales to the same foreground level.
}
\usage{
ssn(x.lumi, targetArray = NULL, scaling = TRUE, bgMethod=c('density', 'mean', 'median', 'none'), fgMethod=c('mean', 'density', 'median'), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.lumi}{ an ExpressionSet inherited object or a data matrix with columns as samples and rows as genes }
  \item{targetArray}{ A target chip is the model for other chips to normalize. It can be a column index, a vector or a LumiBatch object with one sample. }
  \item{scaling}{ determine whether do scaling or just background shift }
  \item{bgMethod}{ optional methods of determining the background level }
  \item{fgMethod}{ optional methods of determining the foreground level }
  \item{\dots}{ other parameters used by \code{\link[stats]{density}} function }
}
\details{
	This function basically adjusts the samples to the same background level and then optionally scales to the same foreground level. The adjustment is based on the raw scale data (For the transformed data, it still estimates the parameters in the raw scale by inverse transformation.). 
	
	Comparing with other normalization methods, like quantile and curve-fitting methods, SSN is a more conservative method. The only assumption is that each sample has the same background levels and the same scale (if do scaling). There are three methods ('density', 'mean' and 'median') for background estimation. If bgMethod is 'none', then the background level will be set as 0, i.e., no background adjustment. For the 'density' bgMethod, it estimates the background based on the mode of probe intensities based on the assumption that the background level intensity is the most frequent value across all the probes in the chip. For the foreground level estimation, it also provides three methods ('mean', 'density', 'median'). For the 'density' fgMethod, it assumes the background probe levels are symmetrically distributed. Then we estimate the foreground levels by taking the intensity mean of all other probes except from the background probes. For the 'mean' and 'median' methods (for both bgMethod and fgMethod), it basically estimates the level based on the mean or median of all probes of the sample. If the fgMethod is the same as bgMethod (except 'density' method), no scaling will be performed.   
}
\value{
Return an object with expression values normalized. The class of the return object is the same as the input object x.lumi.
}
\references{ }
\author{ Pan Du, Simon Lin }
\seealso{ \code{\link{lumiN}} }
\examples{

}
\keyword{ methods }
