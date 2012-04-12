\name{seq2id}
\alias{seq2id}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Transfer a nucleotide sequence as a nuID}
\description{
The nuID (nucleotide universal identifier) is uniquely corresponding to probe sequence. The nuID is also self-identification and error checking
}

\usage{
seq2id(seq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seq}{ a nucleotide sequence composed of A, C, G, T (U). }
}
\details{
  The nuID is a exact mapping of nucleotide sequence based on Base64 encoding scheme. 
A character set A-Z, a-z, 0-9, "\_" and "." is used to represent to the base-64 numbers of 0-63.
The first character of nuID is a checking code, which provide information of both the number of 
padded "A"s at the nucleotide sequence and error checking.
Please refer to reference for more details.
}
\value{
A string represents nuID
}
\references{
Du, P., Kibbe, W.A. and Lin, S.M., "nuID: A universal naming schema of oligonucleotides for Illumina, Affymetrix, and other microarrays", Biology Direct 2007, 2:16 (31May2007).
}
\author{ Pan Du }
\seealso{\code{\link{id2seq}}}
\examples{
seq <- 'ACGTAAATTTCAGTTTAAAACCCCCCG'
id <- seq2id(seq)
id
id2seq(id)
}
\keyword{methods}
