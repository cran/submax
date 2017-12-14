\name{submax}
\alias{submax}
\concept{Causal inference}
\concept{Simultaneous inference}
\concept{Sensitivity analysis}
\concept{Effect modification}
\concept{Subgroup analyses}
\concept{Observational studies}
\concept{rbounds}
\title{
Effect Modification Using the Submax Method in Observational Studies
}
\description{
Effect modification means that the magnitude or stability of a
treatment effect varies with observed covariates.  When there
is effect modification, causal conclusions may be less sensitive
to unmeasured biases in a subgroup in which the treatment effect
is larger or more stable.  The submax or subgroup maximum method
looks at an overall test and subgroup tests, correcting for
multiple testing using the joint distribution of the tests.
The submax method was proposed by Lee et al. (2017).
}
\usage{
submax(y, cmat, gamma = 1, alternative = "greater", alpha = 0.05, rnd = 2, fast=FALSE)
}
\arguments{
  \item{y}{
In general, y is either a matrix with I rows for I matched sets or a vector
of length I for I matched pairs.  The y values are outcomes, perhaps
after scoring (e.g., Huber-Maritz M-scores) or ranking (e.g., Wilcoxon)
to produce a robust test.  If y contains the outcomes themselves, then
the outcomes are permuted in the manner of a permutational t-test.

More precisely, y contains the scores \eqn{q[gij]} discussed in section 2.1
of Lee et al. (2017).  If every matched set gi is a pair, then y may be
a vector of treated-minus-control pair differences, \eqn{q[gi1]-q[gi2]}.

If there are I matched pairs, then y can be either a vector of length I giving
the I treated-minus-control pair differences, or y can be a 2-column matrix
with treated responses in the first column and control responses in the second
column.

If there are I matched sets and the largest matched set contains J individuals,
then y is an I by J matrix with one row for each matched set.  If matched
set i contains one treated individual and k controls, where k is at least 1 and
at most J-1, then y[i,1] is the treated individual's response, y[i,2],...,y[i,k+1]
are the responses of the k controls, and y[i,k+2],...,y[i,J] are equal to NA.

Although y can contain NA's, y[i,1] and y[i,2] must not be NA for every i.  That
is, every matched set must have at least one treated subject and one control.
}
  \item{cmat}{
A matrix with one row for each matched set in y.  For each column of cmat,
a statistical test is done.  Typically, the first column is (1,...,1) and
refers to a test that uses all of the matched sets.  If the second column
is 1 for matched sets of women and 0 for matched sets of men, then the
second test is restricted to matched sets of women.
}
  \item{gamma}{
gamma is the sensitivity parameter \eqn{\Gamma}, where \eqn{\Gamma \ge 1}.  Setting
\eqn{\Gamma = 1} is equivalent to assuming ignorable treatment assignment given
the matched sets, and it performs a within-set randomization test.
}
  \item{alternative}{
If alternative="greater", the null hypothesis of no treatment effect
is tested against the alternative of a treatment effect larger than zero.  If
alternative="less", the null hypothesis of no treatment effect
is tested against the alternative of a treatment effect smaller than 0.
In particular, alternative="less" is equivalent to: (i) alternative="greater",
(ii) y replaced by -y.  See the note for
discussion of two-sided sensitivity analyses.
}
  \item{alpha}{
The global null hypothesis of no effect is tested at
simultaneous level \eqn{\alpha}
in the presence of a bias of at most \eqn{\Gamma}.
}
  \item{rnd}{
The correlation matrix of the dim(cmat)[2] test statistics is
returned, rounded to rnd digits.  The critical.constant is
also rounded to rnd digits.
}
  \item{fast}{
Determines the speed and accuracy of the determination of
the critical.constant used in testing.  fast=TRUE is
faster but less precise, less stable.  fast=FALSE is
slower but more precise, more stable.  See details.
}
}
\details{
The submax procedure is developed by Lee et al. (2017),
and the example reproduces analyses from that paper.

The submax() function rejects the null hypothesis
at level \eqn{\alpha} in the presence of a bias in treatment
assignment of at most \eqn{\Gamma} if maxdeviate is
greater than or equal to critical.constant.

The global null hypothesis of no effect is tested at
simultaneous level \eqn{\alpha}
in the presence of a bias of at most \eqn{\Gamma}.  If the global
null is true, and the bias in treatment assignment is at most
\eqn{\Gamma}, then the probability of falsely rejecting the
global null hypothesis is at most \eqn{\alpha}.
The test looks at the largest of dim(cmat)[2] standardized test
statistics and corrects for multiple testing using their joint
distribution.  The joint distribution is approximated by
a multivariate Normal distribution, so the entire procedure
is a large sample approximation.

The function score() in this package may be helpful in
creating the matrics y and cmat that are arguments of
the submax function.

The sensitivity bound is based on the separable approximation
described in Gastwirth, Krieger and Rosenbaum (2000); see also
Rosenbaum (2007).

The critical.constant is determined rather precisely by a call
to the qmvnorm() function in the mvtnorm package.  The qmvnorm()
function uses random numbers, but this should produce
negligible variability in the critical.constant.  However, if you
call submax() twice, there will be a tiny change in the
critical.constant.  There is
a trade-off between precision and speed, and you can alter
that trade-off -- make submax faster or more precise -- by
editing the call to qmvnorm() in the rcode for submax().
For almost all users, there will be no need to alter the code.
}

\value{
\item{maxdeviate }{The submax or subgroup maximum statistic.  It
is the maximum standardized deviate for the several columns of
cmat.  In Lee et al. (2017), this is \eqn{D[\Gamma max]}}.
\item{critical.constant }{If maxdeviate >= critical.constant,
the global null hypothesis of no treatment effect is
rejected at level \eqn{\alpha} in the presence of a bias
of at most \eqn{\Gamma} }
\item{deviates}{There is one standardized deviate for each
column of cmat, and the maximum of these is maxdeviate.
If deviate[j] > critical.constant, then deviate[j] would
lead to rejection of the global null hypothesis.
}
\item{correlation}{The correlation matrix of the
deviates.  By default, it is printed to two digits
for easy viewing. By changing rnd=2, it can be
produced with additional digits, perhaps for use
in further computations.}
\item{detail}{Reminds the user of the value of
\eqn{\alpha} and \eqn{\Gamma}.}
}
\references{
Cox, D. R. (1977).  The role of signficance tests (with Discussion).
Scand. J. Statist. 4, 49-70.

Gastwirth, J. L., Krieger, A. M. and Rosenbaum, P. R. (2000).
Asymptotic separability in sensitivity analysis.
J. Roy. Statist. Soc. B. 62 545-555. <doi:10.1111/1467-9868.00249>

Lee, K., Small, D. S. and Rosenbaum, P. R. (2017). A new, powerful
approach to the study of effect modification in observational
studies. <arXiv:1702.00525>.

Marcus, R., Eric, P. and Gabriel, K. R. (1976). On closed testing
procedures with special reference to ordered analysis of variance.
Biometrika, 63, 655-660.

Rosenbaum, P. R. (2007). Sensitivity analysis for m-estimates, tests
and confidence intervals in matched observational studies.
Biometrics 63 456-64. (See section 4.)
<doi:10.1111/j.1541-0420.2006.00717.x>
}
\author{
Paul R. Rosenbaum
}
\note{
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

A 2-sided, \eqn{\alpha}-level test may be obtained by performing
two one-sided tests, each at level \eqn{\alpha/2}.  This is a
safe approach, but for \eqn{\Gamma > 1} it is slightly
conservative.  Cox (1977) suggests that we view a two-sided
test as two one-sided tests with a Bonferroni correction
for testing two hypotheses.
}

\examples{
# Reproduces parts of Table 2 of Lee et al. (2017).
data(Active)
submax(Active$delta,Active[,1:7],gamma=1.77,alternative="less",fast=TRUE)
amplify(1.77,c(2,3,4))

# Reproduces the closed-testing analysis in
#Section 4 of Lee et al. (2017)
submax(Active$delta,Active[,c(3,5)],gamma=1.4,alternative="less",fast=TRUE)

# See also the examples for the score() function.
}

\keyword{ htest }
\keyword{ robust }
