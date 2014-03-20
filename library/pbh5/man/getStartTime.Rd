\name{getStartTime}
\alias{getStartTime}
\alias{getAdvanceTime}

\title{
  Computes the Start Time of Polymerization
}
\description{
  'getStartTime' computes the start time, in seconds, of polymerization
  for each base in specified indexes. This dataset is optional so this
  function might not succeed.
}
\usage{
  getStartTime(cmpH5, idx, unit = c("seconds", "frames"))
}
\arguments{
  \item{cmpH5}{
    An object of class \code{PacBioCmpH5}.
  }
  \item{idx}{
    The indices of the alignments to return.
  }
  \item{unit}{
    The unit of the returned values
  }
}
\value{
  'getStartTime' returns a list of numeric vectors, where each vector
  corresponds to one alignment, and each value within that vector
  corresponds to the start time of polymerization for a particular nucleotide.
}
\examples{
  \dontrun{
  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5", package = "pbh5"))
  getStartTime(cmpH5)
  }
}
\keyword{datasets}

