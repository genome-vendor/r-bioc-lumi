\name{example.lumi}
\alias{example.lumi}
\docType{data}
\title{Example LumiBatch object includes example data}
\description{
  Example data as a LumiBatch object which is a subset of Barnes data (Barnes, 2005)
}
\usage{data(example.lumi)}
\format{
A 'LumiBatch' object
}

\details{
The data is from (Barnes, 2005). It used Sentrix HumanRef-8 Expression BeadChip. Two samples "100US" and "95US:5P" (each has two technique replicates) were selected. In order to save space, 8000 genes were randomly selected. As a result, the example data includes 8000 genes, each has 4 measurements. The full data set was included in the Bioconductor Experiment data package lumiBarnes.

The entire data set has been built as a lumiBarnes data object and can be downloaded from Bioconductor Experiment Data.
}
\references{
Barnes, M., Freudenberg, J., Thompson, S., Aronow, B. and Pavlidis, P. (2005) Ex-perimental comparison and cross-validation of the Affymetrix and Illumina gene expression analysis platforms, Nucleic Acids Res, 33, 5914-5923.

The detailed data information can be found at: http://www.bioinformatics.ubc.ca/pavlidis/lab/platformCompare/
}

\examples{
## load the data
data(example.lumi)

## summary of the data
example.lumi

}
\keyword{datasets}
