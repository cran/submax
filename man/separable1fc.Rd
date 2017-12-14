\name{separable1fc}
\alias{separable1fc}
\title{
Computes the Separable Approximation.
}
\description{
Of limited interest to most users, separable1fc() is called by
the main function, submax().
}
\usage{
separable1fc(ymat, gamma = 1)
}
\arguments{
  \item{ymat}{
A matrix of scores produced by mscoref.
}
  \item{gamma}{
The sensitivity parameter \eqn{\Gamma \ge 1}.
}
}
\details{
See Gastwirth, Krieger and Rosenbaum (2000) and
Rosenbaum (2007, section 4) for discussion of the
separable approximation.
}
\value{
\item{tstat }{Vector of length I = dim(ymat)[1] giving the values
of the test statistic in the I matched sets.}
\item{expect }{Vector of length I giving the maximum
expectations in the I matched sets.}
\item{vari }{Vector of length I giving the maximum
variances at the maximum expectations in the I matched sets.}
}

\author{
Paul R. Rosenbaum
}
\note{
This function is similar to the separable1f() function in the
sensitivityfull package.  Unlike that function, separable1fc()
returns the I components for the I matched sets, rather than
computing a summary statistic from them.
}

\references{
Gastwirth, J. L., Krieger, A. M. and Rosenbaum, P. R. (2000).
Asymptotic separability in sensitivity analysis.
J. Roy. Statist. Soc. B. 62 545-555. <doi:10.1111/1467-9868.00249>

Rosenbaum, P. R. (2007). Sensitivity analysis for m-estimates, tests
and confidence intervals in matched observational studies.
Biometrics 63 456-64. (See section 4.)
<doi:10.1111/j.1541-0420.2006.00717.x>
}


\examples{
# The following artificial example computes mscores for a
# full matching, then applies separable1fc() to
# perform a sensitivity analysis.  Compare with
# the example below from the sensitivityfull package.

# The artificial example that follows has I=9
# matched sets.  The first 3 sets have one treated
# individual and two controls with treated subjects
# in column 1.  The next 3 sets are
# matched pairs, with treated subjects in column 1.
# The next 3 sets have one control and two treated
# subjects, with the control in column 1.  Simulated
# from a Normal distribution with an additive effect
# of tau=1.

y<-c(2.2, 1.4, 1.6, 2.4, 0.7, 1.3, 1.2, 0.6, 0.3,
0.5, -0.1, -1.3, -0.3, 0.1, 0.4, 3.0, 1.1, 1.4, -0.8,
0.1, 0.8, NA, NA, NA, 1.1, 0.5, 1.8)
y<-matrix(y,9,3)
treated1<-c(rep(TRUE,6),rep(FALSE,3))

s<-separable1fc(sensitivityfull::mscoref(y,treated1),gamma=2)
1-pnorm((sum(s$tstat)-sum(s$expect))/sqrt(sum(s$vari)))
sensitivityfull::senfm(y,treated1,gamma=2)
s
}
\keyword{ htest }
\keyword{ robust }
