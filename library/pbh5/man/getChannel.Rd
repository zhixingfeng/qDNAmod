\name{getChannel}
\alias{getChannel}
\alias{getChi2}
\alias{getIsPulse}   
\alias{getMaxSignal}  
\alias{getMeanSignal}   
\alias{getMidSignal}  
\alias{getStartFrame}
\alias{getWidthInFrames}
\alias{getMidStdDevSignal}
\alias{getClassifierQV}
\title{
  Pulsecall-level Metrics 
}
\description{
  Retrieve matrices (or vectors) of features stored in the
  pls.h5 files. This information is not necessarily stored in the pls.h5
  file and code should check before depending on its
  existence. Additionally, some of these fields are available in both
  pls and bas files and caution should be taken in ensuring that the
  appropriate value is retrieved. 
}
\usage{
  getChannel(plsH5, holeNumbers = getHoleNumbers(plsH5), ...)
}
\arguments{
  \item{plsH5}{
    An object created by \code{PacBioPlsH5}
  }
  \item{holeNumbers}{
    Which ZMWs to grab.
  }
  \item{...}{
    currently ignored
  }
}
\value{
  A list of numeric vectors or matrices depending on the size of the
  dataset to be retrieved.
}
\details{
  In order to determine the meaning of the particular features consult
  the PacBio bas.h5 documentation.
}
\keyword{datasets}
