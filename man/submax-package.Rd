\name{submax-package}
\alias{submax-package}
\docType{package}
\title{
\packageTitle{submax}
}
\description{
\packageDescription{submax}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{submax}
\packageIndices{submax}
The main function is submax().  Also helpful is score().  See their documentation.
}
\author{
\packageAuthor{submax}

Maintainer: \packageMaintainer{submax}
}
\references{
Lee, K., Small, D. S., & Rosenbaum, P. R. (2017). A new, powerful
approach to the study of effect modification in observational
studies. arXiv preprint arXiv:1702.00525.
}

\keyword{ package }

\examples{
#Reproduces parts of Table 2 of Lee et al. (2017)
data(Active)
submax(Active$delta,Active[,1:7],gamma=1.70,alternative="less")
}
