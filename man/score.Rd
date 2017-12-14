\name{score}
\alias{score}
\concept{M-estimate}
\concept{M-test}
\title{
Creates M-scores for Use by submax().
}
\description{
The score() function is an optional aid in using
the submax() function.  The submax() function may be
used on its own, but the user must then create
the y matrix and the cmat matrix.  The score()
function can be used to create the y matrix and
cmat matrix for use in the submax() function.
It creates Huber-Maritz M-scores.
}
\usage{
score(y, z, mset, x, expandx = FALSE, scale = "closed", inner = 0, trim = 3,
lambda = 1/2, xnames=NULL)
}
\arguments{
  \item{y}{ A vector of responses with no missing data.}
  \item{z}{ Treatment indicator, z=1 for treated, z=0 for control with
  length(z)==length(y).
}
  \item{mset}{ Matched set indicator, 1, 2, ..., sum(z) with
  length(mset)==length(y).  Matched set indicators should be either integers or a factor.
}
  \item{x}{
A matrix of binary covariates.  These are the potential effect modifiers.
If x is a binary vector for one covariate, it will be reshaped into a matrix.
An error will result if x does not have 1 or 0 coordinates.  The matrix
x should length(y) rows.
}
  \item{expandx}{If expandx=FALSE, then x is used as is.  If expandx=TRUE,
then a column of 1's is added to x, as well as 1-x[,j] for each column j.
For example, if x is a single column, 1 for female, 0 for male, then
with expandx=TRUE, x will be replaced by 3 columns, all 1's, female,
and 1-female=male. }
  \item{scale}{scale determines how the observations are scaled
in computing M-scores.  scale must equal "closed"" or "global" or
"interaction".  If you use the mean (or equivalently the total) as
your test statistic (by setting inner=0 and trim=Inf), then scaling
is not needed and the scale parameter is ignored.
If scale=global, then all observations are used in
computing a single scale factor.  This is appropriate when testing
Fisher's global hypothesis of no treatment effect at all.  If
scale=interaction, then the scale is determined separately in
each of the unique groups formed from the interaction of all
of the columns of x.  This can be reasonable if the interaction
groups are not too small.  If x had a column for men, a column for
women, a column for people over 50, a column for people under
50, then there would be 4 interaction groups, such as
men under 50.  If scale=closed, then a single scale factor
is computed using every y[i] such that x[i,j]=1 for at least
one j.  With scale=closed, the score() function will return matrices y and cmat
that only contain the rows i such that x[i,j]=1 for at least
one j.  For instance, if your hypotheses are confined to women,
then the men will be excluded.  score=closed is useful in closed testing,
as described in
Lee et al. (2017).  If x contains a column of 1's, then
scale=closed and scale=global are equivalent.  Of course, x will
contain a column of 1's if you set expandx=TRUE.  You can use
scale=interaction only if sets are exactly matched for
all effect modifiers in x, but this is not required for
scale=global or scale=closed.  If you set scale=interaction
when sets are not exactly matched, you will receive a warning
and scale will be reset to the default of closed.}
  \item{inner}{
inner and trim together define the \eqn{\psi}-function for the M-statistic.
The default values yield a version of Huber's \eqn{\psi}-function, while
setting inner = 0 and trim = Inf uses the mean or total within each matched set.
The \eqn{\psi}-function is an odd function, so \eqn{\psi(w) = -\psi(-w)}.  For
\eqn{w \ge 0}, the \eqn{\psi}-function is \eqn{\psi(w)=0} for \eqn{0 \le w \le } inner,
is \eqn{\psi(w)= } trim for \eqn{w \ge } trim, and rises linearly from 0 to trim
for inner < w < trim.

If uncertain about inner, trim and lambda, then use the defaults.

An error will result unless \eqn{0 \le } inner \eqn{ \le } trim.
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
  \item{xnames}{
If xnames=NULL and x is a matrix or data.frame with column names, then
those names are used to label output.  If xnames is
not null, and x is a matrix, then xnames must be a vector of dim(x)[2]
names, and these names are used to label output.  If x is a vector,
not a matrix, then the output will be easier to read if you give it
a name, xnames="a.name".  This is particularly true if x is a
vector and expandx=TRUE.
}
}
\details{
Taking trim < Inf limits the influence of outliers; see Huber (1981).

Taking trim < Inf and inner = 0 uses Huber's psi function.

Taking trim = Inf and inner = 0 does no trimming and is similar to a mean
or a weighted mean.

Taking inner > 0 often increases design sensitivity; see Rosenbaum (2013).  This
is most evident with matched pairs, where inner=0.5 may be a good choice.

An M-statistic
similar to a trimmed mean is obtained by setting inner=0, trim=1, and
1-lambda to the total amount to be trimmed from both tails.  For example,
inner=0, trim=1, lambda=0.9 trims 10 percent, perhaps 5 percent from
each tail.  Arguably, inner=0, trim=1, and lambda=0.99 is very much
like a mean, but also safer than a mean.

In general, each call to submax() tests a global null hypothesis of
no treatment effect, but does this by looking in several subgroups
defined by pretreatment covariates, that is, by potential effect
modifiers.  We often want to say something about the subgroups
themselves, not the global null hypothesis.  This is possible
by using closed testing (Marcus et al. 1976), which may entail
several calls to
submax().  Before using closed testing, it is suggested that
you read the section in Lee et al. (2017) discussing closed
testing.

Matched sets of unequal size are weighted with a view to
efficiency.  See the documentation for mscorev() and the
references mentioned there.
}
\value{
\item{y }{A matrix of M-scores suitable for use as the y matrix
in the submax() function.}
\item{cmat }{A matrix of comparisons suitable for use as the
cmat matrix in the submax funtion.}
\item{detail }{A data.frame reminding you of the settings that
produced the M-scores.  It contains inner, trim, lambda, scale,
permutational.t, and anyinexact, where permutational.t=TRUE
if you set inner=0 and trim=Inf, and anyinexact=TRUE if some
of the potential effect modifiers are not exactly matched.}
}
\references{

Huber, P. (1981) Robust Statistics.  New York: John Wiley.
(M-estimates based on M-statistics.)

Lee, K., Small, D. S., & Rosenbaum, P. R. (2017). A new, powerful
approach to the study of effect modification in observational
studies. arXiv:1702.00525.

Marcus, R., Eric, P., & Gabriel, K. R. (1976). On closed testing
procedures with special reference to ordered analysis of variance.
Biometrika, 63, 655-660.

Maritz, J. S. (1979). A note on exact robust confidence
intervals for location. Biometrika 66 163--166. (Introduces
exact permutation tests based on M-statistics by redefining
the scaling parameter.)

Rosenbaum, P. R. (2007). Sensitivity analysis for m-estimates, tests
and confidence intervals in matched observational studies.
Biometrics 63 456-64. (R package sensitivitymv) <doi:10.1111/j.1541-0420.2006.00717.x>

Rosenbaum, P. R. (2013). Impact of multiple matched controls on
design sensitivity in observational studies. Biometrics 69 118-127.
(Introduces inner trimming.) <doi:10.1111/j.1541-0420.2012.01821.x>

Rosenbaum, P. R. (2014). Weighted M-statistics with superior design
sensitivity in matched observational studies with multiple controls.
J. Am. Statist. Assoc. 109 1145-1158. (R package sensitivitymw)
<doi:10.1080/01621459.2013.879261>

Rosenbaum, P. R. (2015).  Two R packages for sensitivity analysis
in observational studies.  Observational Studies, v. 1. (Free on-line.)

}
\author{
Paul R. Rosenbaum.
}
\note{
Several other packages use M-scores, so their examples and
documentation may be helpful.  These include sensitivitymv,
sensitivitymw, sensitivitymult, sensitivityfull, and
senstrat.  In particular, sensitivitymw and sensitivitymult
compute sensitivity analyses for confidence intervals.

If you have inexactly matched sets, but you want to
use scale=interaction, then you must discard the
inexactly matched sets before using score().
}


\examples{
data(mercury)
attach(mercury)
# The mercury data has two binary covariates, black and female,
# that will be considered as potential effect modifiers.
# Both black and female are not exactly matched.  Of 397
# matched sets, 72 contain three blacks, 319 contain no
# blacks, 3 contain one black, and 3 contain 2 blacks.
table(table(mset,black)[,2])
# A similar situation arises with females.
table(table(mset,female)[,2])
# When considering females as an effect modifier, only
# sets exactly matched for female are used, etc.  A
# set that is inexact for black may be used when looking
# at females, providing that set is exactly matched for female.


male<-1-female
nonblack<-1-black
everyone<-rep(1,dim(mercury)[1])
x<-cbind(everyone,female,male,black,nonblack)
sc<-score(methylmercury,fish,mset,x)

# At gamma=4, the global null of no effect is
# rejected at alpha=0.05 by every subgroup test
submax(sc$y,sc$cmat,gamma=4,fast=TRUE)

# What does expandx do?
sc<-score(methylmercury,fish,mset,cbind(female,black))
head(sc$cmat)
sc<-score(methylmercury,fish,mset,cbind(female,black),expandx=TRUE)
head(sc$cmat)

# Using exandx with a vector: remember to give it a name.
sc<-score(methylmercury,fish,mset,female,xnames="Female",expandx=TRUE)
head(sc$cmat)

 \dontrun{
# For closed testing, the process is repeated with fewer columns.
# In general, if cmat has L columns, closed testing may require
# up to (2^L)-1 tests.  Here are two of those tests.
  sc<-score(methylmercury,fish,mset,cbind(female,black))
  submax(sc$y,sc$cmat,gamma=4)
# Note that the critical.constant has become smaller, making it
# easier to reject a component hypothesis when fewer hypotheses
# are tested.
  sc<-score(methylmercury,fish,mset,female,xnames="Female")
  submax(sc$y,sc$cmat,gamma=4,fast=TRUE)
# Use of closed testing is discussed in Lee et al. (2017).

# For a two-sided test, change alpha and do 2 tests.
  submax(sc$y,sc$cmat,gamma=4,alpha=0.025,alternative = "greater")
  submax(sc$y,sc$cmat,gamma=4,alpha=0.025,alternative = "less")
# So we reject in the positive direction in all 5 component tests.
  }
detach(mercury)
}

\keyword{ htest }
\keyword{ robust }
