\name{nalevs}
\alias{nalevs}
\title{
Recode variables to 0-1 scale
}
\description{
\code{nalevs} takes as an input any vector and recodes it to range from 0 to 1, to treat specified levels as missing, to treat specified levels as 0, 1, .5, or the mean (weighted or unweighted) of the levels present after coding.
}
\usage{
nalevs(x, naset=NULL, setmid=NULL, set1=NULL, set0=NULL, setmean=NULL, weight=NULL)
}
\arguments{
  \item{x}{A vector to be recoded to range from 0 to 1.}
  \item{naset}{A vector of values of \code{x} to be coded as \code{NA}.}
  \item{setmid}{A vector of values of \code{x} to be recoded to .5.}
  \item{set1}{A vector of values of \code{x} to be recoded to 1.}
  \item{set0}{A vector of values of \code{x} to be recoded to 0.}
  \item{setmean}{A vector of values of \code{x} to be recoded to the mean (if no weight is specified) or weighted mean (if a weight is specified) of values of \code{x} after all recoding.}
  \item{weight}{A vector of weights for \code{x} if weighted means are desired for values listed for \code{setmean}.}
}
\value{
A vector of length equal to that of \code{x} of class \code{numeric}.
}
\author{
Josh Pasek, Assistant Professor of Communication Studies at the University of Michigan (www.joshpasek.com).
}
\examples{
data(anes04)
summary(anes04$age)
summary(nalevs(anes04$age))
table(anes04$educcats)
table(nalevs(anes04$educcats, naset=c(2, 4)))
}
