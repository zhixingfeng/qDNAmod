\name{getZScore}
\alias{getZScore}
\title{
  Compute the Z-Score
}
\description{
  'getZScore' computes the Z-score for each alignment. Z-score is
  generally an internal metric within Pacific Biosciences and might very
  well not exist in external cmp.h5 files.
}
\usage{
getZScore(cmpH5, idx = seq.int(1, nrow(cmpH5)))
}
\arguments{
  \item{cmpH5}{
    An object of class \code{PacBioCmpH5}.
  }
  \item{idx}{
    The indices of the alignments to retrieve.
  }
}
\value{
  'getZScore' is returned as a vector of numbers. 
}
\examples{
  \dontrun{
  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5",
  package = "pbh5"))
  zscore <- getZScore(cmpH5)
  plot(zscore)
  }
}
\keyword{datasets}

