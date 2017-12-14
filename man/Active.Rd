\name{Active}
\alias{Active}
\concept{Observational study}
\concept{Matched pairs}
\concept{Paired survival data}
\docType{data}
\title{
Physical Activity and Survival in NHANES
}
\description{
Physical activity and survival in the NHANES I Epidemiologic Follow-up Study or NHEFS.
This is the example in Lee et al. (2017).  It is patterned after a study by
Davis et al. (1994).  The NHEFS combined the NHANES I study with follow-up for survival.
There are 470 matched pairs consisting of a treated group who were quite inactive
at the time of NHANES I and a matched control group who were very active.
}
\usage{data("Active")}
\format{
  A data frame with 470 observations on the following 12 variables.
  \describe{
    \item{\code{All}}{470 ones, one for each pair}
    \item{\code{male}}{1 for male, 0 for female}
    \item{\code{female}}{1 for female, 0 for male}
    \item{\code{poor}}{1 for income less than 2x poverty level, 0 otherwise}
    \item{\code{notpoor}}{1 for income greater than 2x poverty level, 0 otherwise}
    \item{\code{smoker}}{1 for current smoker, 0 otherwise}
    \item{\code{nonsmoker}}{1 for current nonsmoker, 0 otherwise}
    \item{\code{delta}}{O Brien and Fleming scores for censored matched pairs.  See details.}
    \item{\code{treated.followup.time}}{Death or censoring time for inactive individual}
    \item{\code{treated.censored.time}}{Equals Inf if the inactive individual was censored}
    \item{\code{control.followup.time}}{Death or censoring time for active individual}
    \item{\code{control.censored.time}}{Equals Inf if the active individual was censored}
  }
}
\details{
Pairs were exactly matched for male/female, poor/notpoor and smoker/nonsmoker.
Additionally, pairs were matched for age, white/nonwhite, years of education,
employed or not during the previous 3 months, marital status, alcohol consumption
and dietary quality.  The covariates that were not exactly matched were matched
by minimizing the total Mahalanobis distance within matched pairs.  Table 1 of
Lee et al. (2017) shows covariate balance before and after matching.  These
matching techniques are described in Chapter 8 and Section 9.2 of Rosenbaum (2010).
The matching is the usual kind in epidemiology and biostatistics, that is,
so-called without-replacement matching, in
which no person appears twice.

The example reproduces a row from Table 2 of Lee et al. (2017).

The values in delta in Active are the Prentice-Wilcoxon scores for censored
paired data proposed by O Brien and Fleming (1987).  Specifically, delta
is the \eqn{\Delta} in section 2 of of O Brien and Fleming (1987).
Following their suggestion at the end of their section 2,
data are considered censored at the
earlier censoring time in a matched pair.  These deltas are computed
separately in 8 = 2x2x2 subgroups defined by
male/female x poor/notpoor x smoker/nonsmoker.  Separate
computation of \eqn{\Delta} in subgroups is not needed for the global test of
no treatment effect at all in section 3.2 of Lee et al. (2017),
but it is an aspect of the simultaneous inference by closed
testing in section 4 of Lee et al. (2017).

The NHEFS was, essentially, the NHANES I snapshot survey combined with
follow-up for survival.  Data on mortality and time of death were collected
in four follow-up surveys in 1982-1984, 1986, 1987 and 1992. Tracing of
subjects which enable determination of whether the subject was alive or had
died was high. Ninety six percent of the study population had been successfully
traced at some point through the 1992 follow-up. Tracing rates for each follow-up
ranged from 90 to 94 percent. See Cox et al. (1997).
}
\source{
The data set was constructed by Kwonsang Lee from the NHEFS; see Lee et al. (2017).
The original NHEFS data are publicly available at the NHANES web-page at CDC.
}
\references{
Cox, C. S., Mussolino, M. E., Rothwell, S. T., Lane, M. A., Golden, C. D., Madans, J. H.,
and Feldman, J. J. (1997). Plan and operation of the NHANES I
Epidemiologic Followup Study, 1992. Vital and health
statistics. Ser. 1, Programs and collection procedures, (35), 1-231.

Davis, M. A., Neuhaus, J. M., Moritz, D. J., Lein, D., Barclay, J. D., and
Murphy, S. P. (1994). Health behaviors and survival among middle aged
and older men and women in the NHANES I Epidemiologic Follow-Up Study.
Preventive Medicine, 23, 369-376.

Lee, K., Small, D. S., & Rosenbaum, P. R. (2017). A new, powerful
approach to the study of effect modification in
observational studies. <arXiv:1702.00525>.

O Brien, P. C. and Fleming, T. R. (1987).  A paired Prentice-Wilcoxon
test for censored paired data.  Biometrics, 43, 169-180.
The variable delta in Active is the delta in section 2 of this paper.
Following their suggestion at the end of their section 2,
data are considered censored at the
earlier censoring time in a matched pair.

Rosenbaum, P. R. (2010). Design of Observational Studies.
New York: Springer.
}
\examples{
# The example is from Lee et al. (2017).
data(Active)
submax(Active$delta,Active[,1:7],gamma=1,alternative="less")
}
\keyword{datasets}
