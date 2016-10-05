SampleSizeDiscSurv <- function(power=0.9,alpha=0.025,alternative=c("less","greater"),beta0=0,
                               h0,h1,p0,p1,ties.method=c("efron","breslow","PrenticeGloeckler"),
                               method=c("asymptotic","simulation"),tol,AMV=NULL,nsim=10000,Nvec=NULL,test=c("Wald","Score"))
{
  if (substr(alternative,1,1)=="l") alternative="less" else alternative="greater"
  if (substr(ties.method,1,1)=="P") ties.method="PrenticeGloeckler" else if (substr(ties.method,1,1)=="b")
    ties.method="breslow" else  ties.method="efron"
  if (substr(method,1,1)=="s") method="simulation" else method="asymptotic"
  if (substr(test,1,1)=="W") test="Wald" else test="Score"

  if (abs(beta0)>0 & method=="simulation" & ties.method=="PrenticeGloeckler")
    stop("For simulation method with PrenticeGloeckler ties method, beta0=0 is the only allowed option",
                                                                                  call. = FALSE)

  if (substr(method,1,1)=="a") {
    EZobj=qnorm(power)+qnorm(1-alpha)
    if (alternative=="less") EZobj=-EZobj
    if (is.null(AMV)) AMV=AsympDiscSurv(h0,h1,p0,p1,method=method,tol=1E-12)
    slope.est=(AMV$estimate-beta0)/sqrt(AMV$varn)
    sqrtN=EZobj/slope.est
    if (sqrtN<2) N=0 else N=ceiling(sqrtN^2)
    SSDSl=list(N=N,alternative=alternative,beta0=beta0,ties.method=ties.method,method=method,AMV=AMV,EZobj=EZobj,Nvec=NULL,
               EZvec=NULL,VZvec=NULL,int.est=0,slope.est=slope.est,nsim,test)} else {
    if (is.null(Nvec)) {

      SSDSa=SampleSizeDiscSurv(power,alpha,alternative,beta0,h0,h1,p0,p1,ties.method,
                               method="asymptotic",tol,AMV,nsim,Nvec,test)
      AMV=AMV
      Nvec=ceiling(c(max(c(SSDSa$N/2,SSDSa$N-200)),SSDSa$N+200))
    }
    EZvec=Nvec
    VZvec=Nvec
    lh0=length(h0)
    Zi=rep(0,nsim)

    mprob=c(p0[1:(lh0-1)]-p0[2:lh0],h0[lh0]*p0[lh0],(1-h0[lh0])*p0[lh0],p1[1:(lh0-1)]-p1[2:lh0],h1[lh0]*p1[lh0],(1-h1[lh0])*p1[lh0])
    rateevent=c(h0*p0,h1*p1)/mprob[c(1:lh0,(lh0+2):(2*lh0+1))]
    if ((p0[1]+p1[1])<1) mprob=c(mprob,1-p0[1]-p1[1])
    for (j in 1:length(Nvec)) {
      for (i in 1:nsim) {
        tm=as.vector(rmultinom(1,Nvec[j],mprob)[1:(2*lh0+2)])
        time=rep(c(1:(lh0+1),1:(lh0+1)),tm)
        grp=rep(1,length(time))
        grp[1:sum(tm[1:(lh0+1)])]=0
        tme=rbinom(2*lh0,tm[c(1:lh0,(lh0+2):(2*lh0+1))],rateevent)
        event=rep(0,length(time))
        if (tme[1]>0) event[1:tme[1]]=1
        cstm=cumsum(tm)
        for (k in 2:lh0) if (tme[k]>0) event[(cstm[k-1]+1):(cstm[k-1]+tme[k])]=1
        if (tm[lh0+1]>0) event[(cstm[lh0]+1):(cstm[lh0]+tm[lh0+1])]=1
        for (k in 2:lh0) if (tme[lh0+k-1]>0) event[(cstm[lh0+k-1]+1):(cstm[lh0+k-1]+tme[lh0+k-1])]=1
        if (tm[2*lh0+2]>0) event[(cstm[2*lh0+2]+1):(cstm[2*lh0+2]+tm[2*lh0+2])]=1
        if (ties.method=="PrenticeGloeckler") {
          PG=PrenticeGloeckler.test(time,event,grp,lh0)
          if (test=="Wald") Zi[i]=PG$wald.test else Zi[i]=sign(PG$coefficient)*sqrt(PG$score.test)} else
        {
          time[event==0]=time[event==0]-0.1
          event[time>lh0]=0
          cx1=coxph(Surv(time[time<(lh0-0.5)],event[time<(lh0-0.5)])~grp[time<(lh0-0.5)],ties=ties.method,init=beta0)
          if (test=="Wald") Zi[i]=sign(cx1$coefficients[1]-beta0)*sqrt(cx1$wald.test) else Zi[i]=sign(cx1$coefficients[1]-beta0)*sqrt(cx1$score)
        }
      }
      EZvec[j]=mean(Zi)
      VZvec[j]=var(Zi)
    }
    sqrtN=sqrt(Nvec)
    lm1=lm(EZvec~sqrtN)
    sigma=sqrt(mean(VZvec))
    EZobj=sigma*qnorm(power)+qnorm(1-alpha)
    if (alternative=="less") EZobj=-EZobj
    int.est=as.numeric(lm1$coef[1])
    slope.est=as.numeric(lm1$coef[2])
    N=ceiling(((EZobj-int.est)/slope.est)^2)
    SSDSl=list(N=N,alternative=alternative,beta0=beta0,ties.method=ties.method,method=method,AMV=AMV,EZobj=EZobj,Nvec=Nvec,
               EZvec=EZvec,VZvec=VZvec,int.est=int.est,slope.est=slope.est,nsim,test)
  }
  class(SSDSl)="SSDS"
  return(SSDSl)
}

print.SSDS=function(x, ...) {
  if (x$N<2) cat("The expected value of the test statistic is in the wrong direction for the alternative specified. Sample size not calculated.\n") else
                   cat(paste("Required sample size = ",x$N,".\n",sep=""))
  if (x$slope>0 & x$method=="simulation") cat(paste("For a total sample size N, the expected value of the test statistic is approximately ",
                                                    signif(x$int.est,4),"+",signif(x$slope.est,4)," * sqrt(N)\n",sep="")) else
  if (x$slope<0 & x$method=="simulation") cat(paste("For a total sample size N, the expected value of the test statistic is approximately ",
                                                    signif(x$int.est,4),signif(x$slope.est,4)," * sqrt(N)\n",sep="")) else
                                         cat(paste("For a total sample size N, the expected value of the test statistic is approximately ",
                                                    signif(x$slope.est,4)," * sqrt(N)\n",sep=""))
  cat(paste("The objective is to make the expected value equal to ",signif(x$EZobj,3),".\n",sep=""))
}
