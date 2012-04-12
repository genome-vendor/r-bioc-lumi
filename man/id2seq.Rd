\name{id2seq}
\alias{id2seq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Transfer a nuID as a nucleotide sequence }
\description{
The nuID (nucleotide universal identifier) is uniquely corresponding to probe sequence. The nuID is also self-identification and error checking
}

\usage{
id2seq(id)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{ a nuID (nucleotide universal identifier) }
}
\details{
A reverse of \code{\link{seq2id}}.
Please refer to reference for more details.
}
\value{
a string of nucleotide sequence
}
\references{
Du, P., Kibbe, W.A. and Lin, S.M., "nuID: A universal naming schema of oligonucleotides for Illumina, Affymetrix, and other microarrays", Biology Direct 2007, 2:16 (31May2007).
}
\author{ Pan Du }
\seealso{\code{\link{seq2id}}}
\examples{
seq <- 'ACGTAAATTTCAGTTTAAAACCCCCCG'
id <- seq2id(seq)
id
id2seq(id)
}
\keyword{methods }
