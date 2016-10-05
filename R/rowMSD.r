rowMSD=function(x)
{ n=ncol(x)
  rm=rowMeans(x)
  return(list(rm=rm,rsd=sqrt(rowSums((x-matrix(rep(rm,n),nrow=nrow(x),ncol=n))^2)/(n - 1))))
}
