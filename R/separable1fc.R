separable1fc <-
function (ymat, gamma = 1) 
{
  n <- dim(ymat)[1]
  m <- dim(ymat)[2]
  mumax <- NULL
  varmax <- NULL
  tstat <- NULL
  for (i in 1:n) {
    if (!is.na(ymat[i, 1])) {
      rk <- sort(as.vector(unlist(ymat[i, ])))
      rk <- sort(rk[!is.na(rk)])
      mi <- length(rk)
      mumaxi <- (-Inf)
      varmaxi <- (-Inf)
      for (ai in 1:(mi - 1)) {
        denom <- ai + gamma * (mi - ai)
        mua <- (sum(rk[1:ai]) + gamma * sum(rk[(ai + 
                                                  1):mi]))/denom
        vara <- ((sum(rk[1:ai]^2) + gamma * sum(rk[(ai + 
                                                      1):mi]^2))/denom) - (mua^2)
        if (mua > mumaxi) {
          mumaxi <- mua
          varmaxi <- vara
        }
        else if (mua == mumaxi) 
          varmaxi <- max(varmaxi, vara)
      }
      mumax <- c(mumax, mumaxi)
      varmax <- c(varmax, varmaxi)
      tstat <- c(tstat, as.vector(ymat[i, 1]))
    }
  }
  #tstat <- as.vector(tstat)
  list(tstat=tstat,expect=mumax,vari=varmax)
}
