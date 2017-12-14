submax <-
function (y, cmat,  gamma=1, alternative = "greater", alpha=0.05, rnd=2, fast=FALSE)
{
  #Check and organize input
  stopifnot((alternative == "greater") | (alternative == "less"))
  stopifnot((fast==TRUE)|(fast==FALSE))
  stopifnot(gamma >= 1)
  if (is.data.frame(y)) y <- as.matrix(y)
  if (is.vector(y)) {
    y <- y[!is.na(y)]
    treat <- y/2
    cont <- (-y/2)
    y <- cbind(treat, cont)
  }
  stopifnot(is.matrix(y))
  stopifnot(all(!is.na(as.vector(y[, 1:2]))))
  if (is.data.frame(cmat)) cmat<-as.matrix(cmat)
  if (is.vector(cmat)) cmat<-matrix(cmat,length(cmat),1)
  stopifnot(is.matrix(cmat))
  stopifnot((dim(y)[1])==(dim(cmat)[1]))
  if (alternative == "less") {
    y <- (-y)
  }

  #Compute test
  s<-separable1fc(y,gamma=gamma)
  int<-cmat
  for (j in 1:(dim(cmat)[2])){
    int[,j]<-int[,j]*s$vari
  }
  cv<-t(cmat)%*%int
  ct<-as.vector(t(cmat)%*%s$tstat)
  cmu<-as.vector(t(cmat)%*%s$expect)
  sig<-as.vector(sqrt(diag(cv)))
  dvec<-(ct-cmu)/sig
  co<-cv/outer(sig,sig,"*")
  if (fast) crit<-mvtnorm::qmvnorm(1-alpha,sigma=co,ptol=0.001,maxiter=700)$quantile
  else crit<-mvtnorm::qmvnorm(1-alpha,sigma=co,ptol=0.0001,maxiter=2000)$quantile
  maxdeviate<-max(dvec)
  if (!is.null(colnames(cmat))) names(dvec)<-colnames(cmat)
  detail<-c(alpha,gamma)
  names(detail)<-c("alpha","gamma")
  list(maxdeviate=maxdeviate,critical.constant=round(crit,rnd),
       deviates=dvec,correlation=round(co,rnd),detail=detail)
}
