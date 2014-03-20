\name{getQualityValue}
\alias{getQualityValue}
\alias{getQualityValue,PacBioCmpH5-method}
\alias{getQualityValue,PacBioPlsH5-method}
\alias{getQualityValue,HasBaseAccessors-method}

\alias{getDeletionQV}
\alias{getDeletionQV,PacBioCmpH5-method}
\alias{getDeletionQV,PacBioBasH5-method}
\alias{getDeletionQV,HasBaseAccessors-method}

\alias{getDeletionTag}
\alias{getDeletionTag,PacBioCmpH5-method}
\alias{getDeletionTag,PacBioBasH5-method}
\alias{getDeletionTag,HasBaseAccessors-method}

\alias{getInsertionQV}
\alias{getInsertionQV,PacBioCmpH5-method}
\alias{getInsertionQV,PacBioBasH5-method}
\alias{getInsertionQV,HasBaseAccessors-method}

\alias{getMergeQV}
\alias{getMergeQV,PacBioCmpH5-method}
\alias{getMergeQV,PacBioBasH5-method}
\alias{getMergeQV,HasBaseAccessors-method}

\alias{getSubstitutionQV}
\alias{getSubstitutionQV,PacBioCmpH5-method}
\alias{getSubstitutionQV,PacBioBasH5-method}
\alias{getSubstitutionQV,HasBaseAccessors-method}

\alias{getSubstitutionTag}
\alias{getSubstitutionTag,PacBioCmpH5-method}
\alias{getSubstitutionTag,PacBioBasH5-method}
\alias{getSubstitutionTag,HasBaseAccessors-method}
\title{
  Return Classifier Quality Values
}
\description{
  Return classifier quality values from a class/file which contains
  them. 
}
\usage{
  getQualityValue(h5Obj, ...)
  getInsertionQV(h5Obj, ...)
  getDeletionQV(h5Obj, ...)
  getMergeQV(h5Obj, ...)
  getSubstitutionQV(h5Obj, ...)

  getDeletionTag(h5Obj, ...)
  getSubstitutionTag(h5Obj, ...)
}
\arguments{
  \item{h5Obj}{
    An object of class \code{\link{PacBioCmpH5-class}} or
    \code{\link{PacBioBasH5-class}} 
  }
  \item{\dots}{
    Extra arguments -- See class documentation for \code{h5Obj}
  }
}
\value{
  A list of vectors of quality values.
}
\examples{
  gc()
  require(pbh5)
  cmpH5 <- PacBioCmpH5(system.file("h5_files",  "aligned_reads.cmp.h5",
  package = "pbh5"))
  values <- getQualityValue(cmpH5)
  class(values)
  length(values)
  length(values[[1]])
  head(values[[1]])
}
\details{
  Quality values are  provided in the standard \emph{phred} scaling,
  \deqn{ QV = -10 \log_{10} p_{err} }
  so that a larger QV values indicates a smaller probability of the
  given error type.

  Deletion QVs correspond to the the probability of a deletion event
  occuring \emph{before} the indicated sequence position, and the
  DeletionTag value will be an uppercase ASCII-encoded DNA base (`A',
  `C', `T', `G') or else `N' if the primary pipeline did not identify
  a "most-likely" deleted base.
  
  Quality values are calculated in the primary analysis pipeline based
  on pulse metrics.
}
\keyword{datasets}
