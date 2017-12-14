\name{mscorev}
\alias{mscorev}

\title{
Computes M-scores for Permuational M-tests.
}
\description{
Of limited interest to most users, function mscorev() computes M-scores.
A similar function function is in the package sensitivitymv.
}
\usage{
mscorev(ymat, inner = 0, trim = 3, lambda = 0.5)
}

\arguments{
  \item{ymat}{
If there are I matched sets and the largest matched set contains J individuals,
then y is an I by J matrix with one row for each matched set.  If matched
set i contains one treated individual and k controls, where k is at least 1 and
at most J-1, then y[i,1] is the treated individual's response, y[i,2],...,y[i,k+1]
are the responses of the k controls, and y[i,k+2],...,y[i,J] are equal to NA.

Although y can contain NA's, y[i,1] and y[i,2] must not be NA for every i.  That
is, every matched set must have at least one treated subject and one control.
}
  \item{inner}{
inner and trim together define the \eqn{\psi}-function for the M-statistic.
The default values yield a version of Huber's \eqn{\psi}-function, while
setting inner = 0 and trim = Inf uses the mean within each matched set.
The \eqn{\psi}-function is an odd function, so \eqn{\psi(w) = -\psi(-w)}.  For
\eqn{w \ge 0}, the \eqn{\psi}-function is \eqn{\psi(w)=0} for \eqn{0 \le w \le } inner,
is \eqn{\psi(w)= } trim for \eqn{w \ge } trim, and rises linearly from 0 to trim
for inner < w < trim.

If uncertain about inner, trim and lambda, then use the defaults.

An error will result unless \eqn{0 \le } inner \eqn{ \le } trim.

Taking trim < Inf limits the influence of outliers; see Huber (1981).
Taking trim < Inf and inner = 0 uses Huber's psi function.
Taking trim = Inf does no trimming and is similar to a weighted mean; see TonT.
Taking inner > 0 often increases design sensitivity; see Rosenbaum (2013).
}
  \item{trim}{
inner and trim together define the \eqn{\psi}-function for the M-statistic.
See inner.
}
  \item{lambda}{
Before applying the \eqn{\psi}-function to treated-minus-control differences,
the differences are scaled by dividing by the lambda quantile of all
within set absolute differences.  Typically, lambda = 1/2 for the median.
The value of lambda has no effect if trim=Inf and inner=0.  See Maritz (1979)
for the paired case and Rosenbaum (2007) for matched sets.

An error will result unless 0 < lambda < 1.
}
}

\value{
Generally, a matrix with the same dimensions as ymat containing the M-scores.
}
\references{
Huber, P. (1981) Robust Statistics.  New York: John Wiley.
(M-estimates based on M-statistics.)

Maritz, J. S. (1979). A note on exact robust confidence
intervals for location. Biometrika 66 163--166. (Introduces
exact permutation tests based on M-statistics by redefining
the scaling parameter.)

Rosenbaum, P. R. (2007) Sensitivity analysis for m-estimates,
tests and confidence intervals in matched observational studies.
Biometrics, 2007, 63, 456-464.  <doi:10.1111/j.1541-0420.2006.00717.x>

Rosenbaum, P. R. (2013). Impact of multiple matched controls on
design sensitivity in observational studies. Biometrics 69 118-127.
(Introduces inner trimming.) <doi:10.1111/j.1541-0420.2012.01821.x>

Rosenbaum, P. R. (2015).  Two R packages for sensitivity analysis
in observational studies.  Observational Studies, v. 1. (Free on-line.)
}
\author{
Paul R. Rosenbaum
}

\note{
The example reproduces Table 3 in Rosenbaum (2007).

Matched sets of unequal size are weighted using weights that would
be efficient in a randomization test under a simple model with
additive set and treatment effects and errors with constant variance;
see Rosenbaum (2007).  Specifically, the total score in set (row) i is
divided by the number ni of individuals in row i,
as in expression (8) in Rosenbaum (2007).
}

\examples{
# The example reproduces Table 3 in Rosenbaum (2007).
data(tbmetaphase)
mscorev(tbmetaphase,trim=1)
}

\keyword{ robust }
