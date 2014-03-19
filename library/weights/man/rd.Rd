\name{rd}
\alias{rd}
\title{
Round Numbers To Text With No Leading Zero
}
\description{
Rounds numbers to text and drops leading zeros in the process.
}
\usage{
rd(x, digits=2, add=TRUE, max=(digits+3))
}
\arguments{
  \item{x}{A vector of values to be rounded (must be numeric).}
  \item{digits}{The number of digits to round to (must be an integer).}
  \item{add}{An optional dichotomous indicator for whether additional digits should be added if no numbers appear in pre-set digit level.}
  \item{max}{Maximum number of digits to be shown if \code{add=TRUE}.}
}
\value{
A vector of length equal to that of \code{x} of class \code{character}.
}
\author{
Josh Pasek, Assistant Professor of Communication Studies at the University of Michigan (www.joshpasek.com).
}
\examples{
rd(seq(0, 1, by=.1))
}
