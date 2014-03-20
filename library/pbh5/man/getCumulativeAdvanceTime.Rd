\name{getCumulativeAdvanceTime}
\alias{getCumulativeAdvanceTime}
\title{
  Compute the CumulativeAdvanceTime
}
\description{
  'getCumulativeAdvanceTime' computes the cumulative advance time for
  each idx. Cumulative advance time is \code{cumsum(pw + ipd)}.
}
\usage{
  getCumulativeAdvanceTime(cmpH5, idx, unit = c("seconds", "frames"))
}
\arguments{
  \item{cmpH5}{
    An object of class \code{PacBioCmpH5}.
  }
  \item{idx}{
    The indices of the alignments to retrieve.
  }
  \item{unit}{
    Whether you want seconds or frames.
  }
}
\details{
 Inter-pulse distances are given in seconds.
}
\value{
  A list of vectors of inter-pulse distances. 
}
\seealso{
  \code{\link{PacBioCmpH5-class}}
}
\examples{
  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5",
                       package = "pbh5"))
  ipd <- getCumulativeAdvanceTime(cmpH5)
  logCumulativeAdvanceTime <- lapply(ipd, function(x) log(x))
  plot(density(logCumulativeAdvanceTime[[1]], na.rm = TRUE),
       main = "CumulativeAdvanceTimes")
}
\keyword{datasets}
