\name{getBaselineSigma}
\alias{getBaselineSigma}
\alias{getBaselineLevel}
\alias{getBaselineSigma,PacBioBasH5-method}
\alias{getBaselineLevel,PacBioBasH5-method}
\alias{getBaselineSigma,MultiPart-method}
\alias{getBaselineLevel,MultiPart-method}
\alias{getBaselineSigma,PacBioPlsH5-method}
\alias{getBaselineLevel,PacBioPlsH5-method}

\title{
  Get Baseline Sigma
}
\description{
  Retrieve the baseline sigma dataset from pls.h5 files.
}
\usage{
  getBaselineSigma(h5Obj, ...)
}
\arguments{
  \item{h5Obj}{
    An object of class \code{PacBioPlsH5}
  }
  \item{...}{
   extra arguments - currently ignored.
  }
}
\details{
  Baseline sigma is a 4-column matrix with nrow = number of ZMWs. Each
  value represents the channel-specific standard deviation.
}
\value{
  A 4-column matrix
}
\seealso{
  \code{\link{PacBioPlsH5-class}}
}
\keyword{datasets}

