\name{probeID2nuID}
\alias{probeID2nuID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mapping Illumina ProbeID as nuID }
\description{
Mapping Illumina ProbeID as nuID
}
\usage{
probeID2nuID(probeID, lib = "lumiHumanV1")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probeID}{ a vector of Illumina ProbeID }
  \item{lib}{ a nuID indexed Illumina annotation library }
}
\details{
}
\value{
return a list with the same length as the input TargetID vector.
}
\references{
Du, P., Kibbe, W.A. and Lin, S.M., "nuID: A universal naming schema of oligonucleotides for Illumina, Affymetrix, and other microarrays", Biology Direct 2007, 2:16 (31May2007).
}
\author{ Pan Du }
\seealso{ See Also  \code{\link{nuID2probeID}} }
\examples{
if (require(lumiHumanV1)) {
    probeID2nuID('0001240020', lib='lumiHumanV1')
}

}
\keyword{ methods }
