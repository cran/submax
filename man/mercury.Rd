\name{mercury}
\alias{mercury}
\docType{data}
\concept{Observational study}
\concept{Matched study}
\title{
NHANES Mercury/Fish Data
}
\description{
Data from NHANES 2009-2010.  397 treated people who ate at least 15 servings of fish or shellfish during the previous month are matched to two controls who ate at most one serving of fish or shellfish.  The values in methylmercury record the level of methylmercury in blood in mu-g/dl.
}
\usage{data("mercury")}
\format{
  A data frame with 1191 observations on the following 6 variables.
  \describe{
    \item{\code{SEQN}}{NHANES 2009-2010 id number}
    \item{\code{methylmercury}}{Methylmercury in blood in mu-g/dl}
    \item{\code{fish}}{1 if ate >= 15 servings of fish or shellfish, 0 if <=1 serving}
    \item{\code{mset}}{Matched set indicator, 1,...,397.}
    \item{\code{female}}{1 if female, 0 if male}
    \item{\code{black}}{1 if black, 0 otherwise}
  }
}
\details{
Sets were matched 2-to-1 for for age, sex, ratio of household income to
the poverty level, education, ethnic group (black, Hispanic, or
other), and cigarette smoking.  A table showing covariate balance
after matching is in Rosenbaum (2014, Table 1).
}
\source{
From NHANES 2009-2010, publicly available at the NHANES web page at CDC.
}
\references{
Rosenbaum, P. R. (2014) Weighted M-statistics with superior design sensitivity
in matched observational studies with multiple controls. Journal of the
American Statistical Association, 2014. <doi:10.1080/01621459.2013.879261>
}
\examples{
data(mercury)
boxplot(mercury$methylmercury~mercury$fish)
}
\keyword{datasets}
