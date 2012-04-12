\name{nuID2probeID}
\alias{nuID2probeID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mapping nuID into Illumina ProbeID }
\description{
  Mapping nuID into Illumina ProbeID
}
\usage{
nuID2probeID(nuID, lib = "lumiHumanV1")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nuID}{ a vector of nuID }
  \item{lib}{ a nuID indexed Illumina annotation library }
}
\details{
}
\value{
return a list with the same length as the input nuID vector.
}
\references{
Du, P., Kibbe, W.A. and Lin, S.M., "nuID: A universal naming schema of oligonucleotides for Illumina, Affymetrix, and other microarrays", Biology Direct 2007, 2:16 (31May2007).
}
\author{ Pan Du }
\seealso{ \code{\link{probeID2nuID}} }
\examples{
if (require(lumiHumanV1)) {
    nuID2probeID("B2J6WGhV.RevOJYff4", lib='lumiHumanV1')
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
