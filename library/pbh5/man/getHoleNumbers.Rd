\name{getHoleNumbers}
\alias{getHoleNumbers}
\alias{getHoleNumbers,PacBioCmpH5-method}
\alias{getHoleNumbers,MultiPart-method}
\alias{getHoleNumbers,PacBioBasH5-method}
\title{
  Get hole numbers Stored in a PacBio HDF5 File.
}
\description{
  Retrieve the ZMW hole numbers stored in the bas.h5 or pls.h5 file.
}
\usage{
  getHoleNumbers(h5Obj, ...)
}
\arguments{
  \item{h5Obj}{
    A class which stores hole numbers, a.k.a., ZMWs (Zero-mode
  waveguides).
}
\item{...}{
  Currently ignored
  }
}
\value{
  A vector of integers representing hole numbers.
}
\keyword{programming}
