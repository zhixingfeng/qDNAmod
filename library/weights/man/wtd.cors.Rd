\name{wtd.cors}
\alias{wtd.cors}
\title{
Produces weighted correlations quickly using C.
}
\description{
\code{wtd.cors} produces a Pearsons correlation coefficient comparing two variables or matrices.  
}
\usage{
wtd.cors(x, y=NULL, weight=NULL)
}
\arguments{
  \item{x}{
\code{x} should be a matrix or vector which the researcher wishes to correlate with \code{y}.
}
  \item{y}{
\code{y} should be a numerical vector or matrix which the researcher wishes to correlate with \code{x}.  If \code{y} is NULL, \code{x} will be used instead
}
  \item{weight}{
\code{weight} is an optional vector of weights to be used to determining the weighted mean and variance for calculation of the correlations.
}
}
\value{
A matrix of the estimated correlation coefficients.
}
\author{
Marcus Schwemmle at GfK programmed the C code, R wrapper by Josh Pasek, Assistant Professor of Communication Studies at the University of Michigan (www.joshpasek.com).
}
\seealso{
\code{\link{wtd.cor}}
\code{\link{stdz}}
\code{\link{wtd.t.test}}
\code{\link{wtd.chi.sq}}
}
\examples{
test <- c(1,1,1,1,1,1,2,2,2,3,3,3,4,4)
t2 <- rev(test)
weight <- c(.5,.5,.5,.5,.5,1,1,1,1,2,2,2,2,2)

wtd.cors(test, t2)
wtd.cors(test, t2, weight)
}
\keyword{ ~correlation }
\keyword{ ~Pearson }
\keyword{ ~weights }
