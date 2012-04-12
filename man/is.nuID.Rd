\name{is.nuID}
\alias{is.nuID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ nuID self-identification }
\description{
  Self-identify nuID (nucleotide universal identifier) by verify the check code value and the checksum value
}
\usage{
is.nuID(id)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{ nuId or other string }
}

\value{
  Return TRUE if id is a nuID, or else return FALSE.
}
\references{
Du, P., Kibbe, W.A. and Lin, S.M., "nuID: A universal naming schema of oligonucleotides for Illumina, Affymetrix, and other microarrays", Biology Direct 2007, 2:16 (31May2007).
}
\author{ Pan Du }

\seealso{  \code{\link{seq2id}}, \code{\link{id2seq}} }
\examples{
## check the function using a random sequence
id <- 'adfasdfafd'
is.nuID(id)		# FALSE

## check the function using a read nuID
seq <- 'ACGTAAATTTCAGTTTAAAACCCCCCG'
id <- seq2id(seq)
is.nuID(id)		# TRUE

}
\keyword{ methods }
