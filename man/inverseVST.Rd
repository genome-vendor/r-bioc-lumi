\name{inverseVST}
\alias{inverseVST}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Inverse VST transform }
\description{
 Inverse transform of VST (variance stabilizing transform), see \code{\link{vst}}.
}
\usage{
inverseVST(x, fun = c('asinh', 'log'), parameter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a VST transformed LumiBatch object or a numeric matrix or vector}
  \item{fun}{ function used in VST transform }
  \item{parameter}{ parameter of VST function }
}
\details{
Recover the raw data from VST transformed data returned by \code{\link{vst}}.
This function can be directly applied to the VST transformed or VST + RSN normalized LumiBatch object to reverse transform the data to the original scale. 
}
\value{
Return the raw data before VST transform
}
\references{
 Lin, S.M., Du, P., Kibbe, W.A.,  "Model-based Variance-stabilizing Transformation for Illumina Mi-croarray Data", submitted
 }
\author{ Pan Du }

\seealso{ \code{\link{vst}} }
\examples{

## load example data
data(example.lumi)

## get the gene expression mean for one chip
u <- exprs(example.lumi)[,1]
## get the gene standard deviation for one chip
std <- se.exprs(example.lumi)[,1]

## do variance stabilizing transform
transformedU <- vst(u, std)

## do inverse transform and recover the raw data
parameter <- attr(transformedU, 'parameter')
transformFun <- attr(transformedU, 'transformFun')
recoveredU <- inverseVST(transformedU, fun=transformFun, parameter=parameter)

## compare with the raw data
print(u[1:5])
print(recoveredU[1:5])

## do inverse transform of the VST + RSN processed data
lumi.N <- lumiExpresso(example.lumi[,1:2])
## Inverse transform. 
## Note: as the normalization is involved, the processed data will be different from the raw data.
lumi.N.raw <- inverseVST(lumi.N)

}
\keyword{ methods }
