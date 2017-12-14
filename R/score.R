score<-function(y, z, mset, x, expandx = FALSE, scale = "closed", inner = 0, trim = 3,
                lambda =1/2, xnames=NULL){
  #Check input
  stopifnot(is.logical(expandx))
  stopifnot((scale=="closed")|(scale=="global")|(scale=="interaction"))
  stopifnot(is.vector(y))
  stopifnot(is.vector(z))
  stopifnot(is.vector(mset))
  if (is.vector(x)) x<-matrix(x,length(x),1)
  if (is.data.frame(x)) x<-as.matrix(x)
  if (!is.null(xnames)){
    stopifnot(length(xnames)==(dim(x)[2]))
    colnames(x)<-xnames
  }
  stopifnot((inner>=0)&(inner<=trim))
  stopifnot((lambda>0)&(lambda<1))
  stopifnot((length(z)==length(y)))
  stopifnot((length(z)==length(mset)))
  stopifnot(length(z)==(dim(x)[1]))
  stopifnot(all(!is.na(y)))
  stopifnot(all((z==0)|(z==1))) #z is 1 for treated, 0 for control
  stopifnot(all((as.vector(x)==1)|(as.vector(x)==0))) # x must be binary
  tbcheck<-table(z,mset)
  ck<-all(tbcheck[2,]==1)&all(tbcheck[1,]>=1)
  if (!ck){
    warning("Every matched set must contain one treated subject and at least one control.")
    stopifnot(ck)
  }

  #Convert y to matrix
  mset<-as.integer(mset)
  o<-order(mset,1-z)
  y<-y[o]
  z<-z[o]
  x<-x[o, , drop=FALSE]
  mset<-mset[o]
  tb<-table(mset)
  nset<-length(tb)
  setsize<-max(tb)

  makeymat<-function(yj){
    ymat<-matrix(NA,nset,setsize)
    m<-0
    for (i in 1:nset){
      ymat[i,1:tb[i]] <- yj[(m+1):(m+tb[i])]
      m<-m+tb[i]
    }
    ymat
  }

  ymat<-makeymat(y)

  # If expandx, add 1-x[,j] and 1 as additional columns
  if (expandx){
    K1<-dim(x)[2]
    for (j in 1:K1){
      x<-cbind(x,1-x[,j])
      colnames(x)[K1+j]<-paste("Not",colnames(x)[j])
    }
    All<-rep(1,length(y))
    x<-cbind(All,x)
  }

  K<-dim(x)[2]

  #Which matched sets are exactly matched?
  exact<-matrix(FALSE,nset,K)
  cmat<-matrix(0,nset,K)
  if (!is.null(colnames(x))) colnames(cmat)<-colnames(x)
  anyinexact<-FALSE

  for (j in 1:K){
    xj<-makeymat(as.vector(x[,j]))
    ex<-round(apply(xj,1,stats::sd,na.rm=TRUE),3)==0
    if (!all(ex)) anyinexact<-TRUE
    exact[,j]<-ex
    cmat[ex,j]<-as.vector(xj[ex,1])
  }

  if (anyinexact) message("Some matched sets are not exactly matched for some effect modifiers.")

  if (anyinexact&(scale=="interaction")) {
    warning("scale=interaction is not a valid option if some matched sets are not exactly matched.
            This option was reset to the default, scale=closed.  Alternatively,
            you could discard any matched set that is not exactly matched for
            all of x, then use scale=interaction.")
    scale<-"closed"
  }

  if ((trim==Inf)&(inner==0)) ys<-ymat
  else if (scale=="global") {
    ys<-mscorev(ymat, inner=inner, trim=trim, lambda=lambda)
  }
  else if (scale=="interaction") {
    # Make the intersection groups g
    ys<-matrix(NA,dim(ymat)[1],dim(ymat)[2])
    g<-factor(as.vector(cmat[,1]))
    for (j in 2:K){
      g<-g:factor(as.vector(cmat[,j]))
    }
    g<-as.integer(g)
    ug<-sort(unique(g))
    tbg<-table(g)
    if (min(tbg)<3) warning("The smallest interaction group g has fewer than 3 matched sets.  Perhaps
                             reconsider scale = interaction.")
    for (gi in g){
      who<-(g==gi)
      ys[who,]<-mscorev(ymat[who, , drop=FALSE], inner=inner, trim=trim, lambda=lambda)
    }
  }
  else if (scale=="closed"){
    used<-apply(cmat,1,max)==1
    ymat<-ymat[used,]
    cmat<-cmat[used, , drop=FALSE]
    ys<-mscorev(ymat,inner=inner,trim=trim,lambda=lambda)
  }
  permutational.t<-(inner==0)&(trim==Inf)
  if (permutational.t) scale<-NA
  detail<-data.frame(inner,trim,lambda,scale,permutational.t,anyinexact)
  list(y=ys,cmat=cmat,detail=detail)
}
