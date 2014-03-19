\name{starmaker}
\alias{starmaker}
\title{
Produce stars from p values for tables.
}
\description{
Recodes p values to stars for use in tables.
}
\usage{
starmaker(x, p.levels=c(.001, .01, .05, .1), symbols=c("***", "**", "*", "+"))
}
\arguments{
  \item{x}{A vector of p values to be turned into stars (must be numeric).}
  \item{p.levels}{A vector of the maximum p value for each symbol used (p<p.level).}
  \item{symbols}{A vector of the symbols to be displayed for each p value.}
}
\value{
A vector of length equal to that of \code{x} of class \code{character}.
}
\author{
Josh Pasek, Assistant Professor of Communication Studies at the University of Michigan (www.joshpasek.com).
}
\examples{
starmaker(seq(0, .15, by=.01))
cbind(p=seq(0, .15, by=.01), star=starmaker(seq(0, .15, by=.01)))
}
