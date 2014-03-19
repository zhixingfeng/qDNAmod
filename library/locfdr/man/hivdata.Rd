\name{hivdata}
\docType{data}
\alias{hivdata}
\title{HIV data set}
\description{
  The data comprises 7680 \eqn{z}-values, each relating to a
  two-sample \eqn{t}-test. The test compares gene expression values for
  4 HIV patients with values for 4 normal subjects; the \eqn{t}-score
  \code{T[i]} for gene \eqn{i} has been transformed to a normal scale,
  \code{z[i] = qnorm(pt(T[i], df=6))}, so that the
  \code{z[i]}'s theoretically would have a
  standard \eqn{N(0,1)} distribution under the null hypothesis. The
  original experiment is described in van't Wout et. al. (2003).
}
\format{A vector containing 7680 \eqn{z}-values}
\usage{data(hivdata)}
\references{van't Wout, et. al., Cellular gene 
  expression upon human immuno-deficiency virus type 1 infection of
  CD4+-T-Cell lines, Journal ofVirology 77, 1392-1402.}
\keyword{datasets}
