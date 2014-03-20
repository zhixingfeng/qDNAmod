\name{restrictCmpH5}
\alias{restrictCmpH5}

\title{
  Limit a cmpH5 to a subset of the reads.
}
\description{
  Construct a new cmpH5 object representing the subset of reads
  specified in mask.
}
\usage{
  restrictCmpH5(cmpH5, mask = rep(TRUE, nrow(cmpH5)))
}

\arguments{
  \item{cmpH5}{An object of type PacBioCmpH5}
  \item{mask}{A vector or booleans.}
}
\details{
  Some of the the finer details of the cmpH5 object are not changed. For
  instance, the sorting facility in cmpH5 supports an index for
  performance, by filtering out certain reads the accuracy of these
  indexes is off, which will make certain operations slower.
}
\value{
  A new cmpH5 file.
}
\examples{
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5", package = "pbh5"))
  nrow(restrictCmpH5(cmpH5, getReadLength(cmpH5) < 500))
}
\keyword{datasets}

