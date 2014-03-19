\name{locfdr}

\alias{locfdr}

\title{Local False Discovery Rate Calculation}

\description{
  Compute local false discovery rates, following the definitions and
  description in references listed below.}

\usage{
locfdr(zz, bre = 120, df = 7, pct = 0, pct0 = 1/4, nulltype = 1, type =
0, plot = 1, mult, mlests, main = " ", sw = 0)
}

\arguments{
  \item{zz}{A vector of summary statistics, one for each case under
  simultaneous consideration.  The
  calculations  assume a large number of cases, say
  \code{length(zz)} exceeding 200.  Results may be improved by transforming
  zz so that its elements are theoretically distributed
  as \eqn{N(0,1)} under the null hypothesis.  See the locfdr vignette
  for tips on creating zz.}

\item{bre}{Number of breaks in the discretization of the \eqn{z}-score axis,
        or a vector of breakpoints fully describing the
        discretization.  If \code{length(zz)} is small, such as when the
        number of cases is less than about 1000, set bre to a number
        lower than the default of 120.}

\item{df}{Degrees of freedom for fitting the estimated
        density \eqn{f(z)}.}
      
\item{pct}{Excluded tail proportions of \eqn{zz}'s when fitting
	   \eqn{f(z)}. \code{pct=0} includes full range of \eqn{zz}'s. pct can
	   also be a 2-vector, describing the fitting range.}

\item{pct0}{Proportion of the \eqn{zz} distribution used in fitting the
	    null density \eqn{f0(z)} by central matching.  If a 2-vector,
	    e.g. \code{pct0=c(0.25,0.60)}, the range [pct0[1],
	    pct0[2]] is used.  If a scalar, [pct0, 1-pct0] is used.}

\item{nulltype}{Type of null hypothesis assumed in estimating \eqn{f0(z)},
	for use in the fdr calculations.  0 is the theoretical null
	\eqn{N(0,1)}, 1 is maximum likelihood estimation, 2 is central
	matching estimation, 3 is a split normal version of 2.}

\item{type}{Type of fitting used for \eqn{f}; 0 is a natural spline, 1 is
	a polynomial, in either case with degrees of freedom df [so
	total degrees of freedom including the intercept is \code{df+1}.]}

\item{plot}{Plots desired.  0 gives no plots. 1 gives single
	    plot showing the histogram of \eqn{zz} and fitted
	    densities \eqn{f} and \eqn{p0*f0}.
	    2 also gives plot of fdr, and the right and
	    left tail area Fdr curves.  3 gives instead the f1 cdf
	    of the estimated fdr curve; plot=4 gives all three plots.}

\item{mult}{Optional scalar multiple (or vector of multiples) of the
	    sample size for calculation of the corresponding
	    hypothetical Efdr value(s).}
	  
\item{mlests}{Optional vector of initial values for (delta0, sigma0) in
  the maximum likelihood iteration.}

\item{main}{Main heading for the histogram plot when \code{plot>0}.}

\item{sw}{Determines the type of output desired.  2 gives a list
  consisting of the last 5 values listed under Value below.
  3 gives the square matrix of dimension bre-1 representing the influence
	  function of log(fdr).  Any other value
	  of sw returns a list consisting of the first 5 (6 if mult is
	  supplied) values listed below.}
}
      

\details{
  See the locfdr vignette for details and tips.
}

\value{
\item{fdr}{the estimated local false discovery rate for each case,
	using the selected type and nulltype.}

\item{fp0}{the estimated parameters delta (mean of f0), sigma
	(standard deviation of f0), and p0, along with
	their standard errors.}

\item{Efdr}{the expected false discovery rate for the non-null cases,
	a measure of the experiment's power as described in Section 3
	of the second reference.  Overall Efdr and right and left values are
	given, both for the specified nulltype and for nulltype 0.  If
	\code{nulltype==0}, values are given for nulltypes 1 and 0.}

\item{cdf1}{a 99x2 matrix giving the estimated cdf of fdr under the
	non-null distribution f1. Large values of the cdf for small fdr
	values indicate good power; see Section 3 of the second
	reference.  Set plot to 3 or 4 to see the cdf1 plot.}

\item{mat}{A matrix of estimates of \eqn{f(x)},
	\eqn{f0(x)}, \eqn{fdr(x)}, etc. at the \eqn{bre-1} midpoints "x"
	of the break discretization, convenient for comparisons and
	plotting.  Details are in the locfdr vignette.}

\item{z.2}{the interval along the zz-axis outside of which $fdr(z)<0.2$,
	the locations of the yellow triangles in the histogram plot.  If
	no elements of zz on the left or right satisfy the criterion,
	the corresponding element of z.2 is NA.}

\item{call}{the function call.}

\item{mult}{If the argument mult was supplied, vector of the ratios of
	    hypothetical Efdr for the supplied multiples of the sample
	    size to Efdr for the actual sample size.}

\item{pds}{The estimates of p0, delta, and sigma.}

\item{x}{The bin midpoints.}

\item{f}{The values of \eqn{f(z)} at the bin midpoints.}

\item{pds.}{The derivative of the estimates of p0, delta, and sigma with
	    respect to the bin counts.}

\item{stdev}{The delta-method estimates of the standard deviations of
  the p0, delta, and sigma estimates.}

}

\references{
  Efron, B. (2004) "Large-scale simultaneous hypothesis testing: the choice of
  a null hypothesis", Jour Amer Stat Assoc, \bold{99}, pp. 96--104
  
  Efron, B. (2006) "Size, Power, and False Discovery Rates"
  
  Efron, B. (2007) "Correlation and Large-Scale Simultaneous Significance
  Testing", Jour Amer Stat Assoc, \bold{102}, pp. 93--103
  
  \url{http://www-stat.stanford.edu/~brad/papers/}
}


\author{Bradley Efron, Brit B. Turnbull, and Balasubramanian Narasimhan}

\examples{
## HIV data example
data(hivdata)
w <- locfdr(hivdata)
}

\keyword{htest}% at least one, from doc/KEYWORDS
\keyword{models}% __ONLY ONE__ keyword per line
	     
