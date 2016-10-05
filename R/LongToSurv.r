rowMSD=function (x)
{n=ncol(x)
rm=rowMeans(x)
return(list(rm=rm,rsd=sqrt(rowSums((x-matrix(rep(rm,n),nrow=nrow(x),ncol=n))^2)/(n - 1))))
}

LongToSurv=function(M,V,L,U,time,p0f,p1f=NULL,method=c("simulation","asymptotic"),
                    conf.type=c("scheduled","unscheduled","none"),nsim=100000) {
  #finds the survival curves for confirmed events defined by a longitudinal model
  #the longitudinal model includes parameters defined by the vector of time points time
  #and the parameter vector theta
  #p0f and p1f are functions that give the vector of probabilities of crossing the threshold at each time point
  #the covariates included in the model are truncated multivariate normal with mean M, covariance V,
  #lower truncation vector L, and upper truncation vector U

  if (substr(conf.type,1,1)=="n") conf.type="none" else if (substr(conf.type,1,1)=="u")
    conf.type="unscheduled" else  conf.type="scheduled"

    if (conf.type=="scheduled") {
      ltm1=length(time)-1
      acc.rate=as.numeric(pmvnorm(L,U,M,sigma=V)) #acceptance rate

      if (substr(method,1,1)=="s") {
        x=rmvnorm(ceiling(nsim/acc.rate),M,sigma=V)
        x=x[colSums(apply(x,1,function(x) x>L & x<U))==length(L),]
        n=nrow(x)
        if (is.null(p1f)) {
          S0fss = function (x,t,p0f,ltm1) {
            #this function calculates the survival curve for a fixed value of the covariate vector x
            p0=p0f(x,t)
            S0=rep(0,ltm1)
            S0[1]=1-p0[1]*p0[2]
            S0[2]=S0[1]-p0[2]*p0[3]+p0[1]*p0[2]*p0[3]
            for (i in 3:ltm1) S0[i]=S0[i-1]*(1-p0[i+1])+S0[i-2]*p0[i+1]*(1-p0[i])
            return(S0)
          }
          S0MSD=rowMSD(apply(x,1,S0fss,t=time,p0f=p0f,ltm1=ltm1))
          LTSlist=list(time=time,S0=c(S0MSD$rm[1:ltm1],S0MSD$rm[ltm1]),
                       S0err=c(S0MSD$rsd[1:ltm1],S0MSD$rsd[ltm1])/sqrt(n),
                       accept.rate=acc.rate)
        } else {
          S01fss = function (x,t,p0f,p1f,ltm1) {
            #this function calculates the survival curves for a fixed value of the covariate x
            p0=p0f(x,t)
            p1=p1f(x,t)
            S0=rep(0,ltm1)
            S1=S0
            S0[1]=1-p0[1]*p0[2]
            S0[2]=S0[1]-p0[2]*p0[3]+p0[1]*p0[2]*p0[3]
            S1[1]=1-p1[1]*p1[2]
            S1[2]=S1[1]-p1[2]*p1[3]+p1[1]*p1[2]*p1[3]
            for (i in 3:ltm1) {
              S0[i]=S0[i-1]*(1-p0[i+1])+S0[i-2]*p0[i+1]*(1-p0[i])
              S1[i]=S1[i-1]*(1-p1[i+1])+S1[i-2]*p1[i+1]*(1-p1[i])}
            return(c(S0,S1))
          }
          S01MSD=rowMSD(apply(x,1,S01fss,t=time,p0f=p0f,p1f=p1f,ltm1=ltm1))
          LTSlist=list(time=time,S0=c(S01MSD$rm[1:ltm1],S01MSD$rm[ltm1]),
                       S0err=c(S01MSD$rsd[1:ltm1],S01MSD$rsd[ltm1])/sqrt(n),
                       S1=c(S01MSD$rm[(ltm1+1):(2*ltm1)],S01MSD$rm[2*ltm1]),
                       S1err=c(S01MSD$rsd[(ltm1+1):(2*ltm1)],S01MSD$rsd[2*ltm1])/sqrt(n),
                       accept.rate=acc.rate)
        }
      } else {
        sdV=sqrt(diag(V))
        #the integration is done after the covariates are transformed to lie in a subset of the unit hypercube by using the probit transformation
        L2=pnorm((L-M)/sdV)
        U2=pnorm((U-M)/sdV)
        if (is.null(p1f)) {
          S0fsa = function (u,t,M,corV,sdV,p0f,ltm1) {
            #this function calculates the survival curves for a fixed value of the coordinatewise probit-transformed covariate vector u.
            v=qnorm(u)
            p0=p0f(v*sdV+M,t)
            S0=rep(0,ltm1)
            S0[1]=1-p0[1]*p0[2]
            S0[2]=S0[1]-p0[2]*p0[3]+p0[1]*p0[2]*p0[3]
            for (i in 3:ltm1) S0[i]=S0[i-1]*(1-p0[i+1])+S0[i-2]*p0[i+1]*(1-p0[i])
            return(S0*dmvnorm(v,sigma=corV)*exp(sum(v^2)/2))
          }
          i1=adaptIntegrate(S0fsa,L2,U2,t=time,M=M,corV=cov2cor(V),sdV=sdV,p0f=p0f,ltm1=ltm1,tol=1e-05,fDim=ltm1,maxEval=100000,absError=0)
          S0=i1$integral*(2*pi)^(length(sdV)/2)/acc.rate
          S0err=i1$error*S0
          LTSlist=list(time=time,S0=c(S0[1:ltm1],S0[ltm1]),S0err=c(S0err[1:ltm1],S0err[ltm1]),
                       accept.rate=acc.rate)
        } else {
          S01fsa = function (u,t,M,corV,sdV,p0f,p1f,ltm1) {
            #this function calculates the survival curves for a fixed value of the coordinatewise probit-transformed covariate vector u.
            v=qnorm(u)
            p0=p0f(v*sdV+M,t)
            p1=p1f(v*sdV+M,t)
            S0=rep(0,ltm1)
            S1=S0
            S0[1]=1-p0[1]*p0[2]
            S0[2]=S0[1]-p0[2]*p0[3]+p0[1]*p0[2]*p0[3]
            S1[1]=1-p1[1]*p1[2]
            S1[2]=S1[1]-p1[2]*p1[3]+p1[1]*p1[2]*p1[3]
            for (i in 3:ltm1) {
              S0[i]=S0[i-1]*(1-p0[i+1])+S0[i-2]*p0[i+1]*(1-p0[i])
              S1[i]=S1[i-1]*(1-p1[i+1])+S1[i-2]*p1[i+1]*(1-p1[i])}
            return(c(S0,S1)*dmvnorm(v,sigma=corV)*exp(sum(v^2)/2))
          }
          i1=adaptIntegrate(S01fsa,L2,U2,t=time,M=M,corV=cov2cor(V),sdV=sdV,p0f=p0f,p1f=p1f,ltm1=ltm1,tol=1e-05,fDim=2*ltm1,maxEval=100000,absError=0)
          S01=i1$integral*(2*pi)^(length(sdV)/2)/acc.rate
          S01err=i1$error*S01
          LTSlist=list(time=time,S0=c(S01[1:ltm1],S01[ltm1]),
                       S0err=c(S01err[1:ltm1],S01err[ltm1]),
                       S1=c(S01[(ltm1+1):(2*ltm1)],S01[2*ltm1]),
                       S1err=c(S01err[(ltm1+1):(2*ltm1)],S01err[2*ltm1]),accept.rate=acc.rate)
        }
      }
    } else
      {
        ltm=length(time)
        acc.rate=as.numeric(pmvnorm(L,U,M,sigma=V)) #acceptance rate
        if (substr(method,1,1)=="s") {
          x=rmvnorm(ceiling(nsim/acc.rate),M,sigma=V)
          x=x[colSums(apply(x,1,function(x) x>L & x<U))==length(L),]
          n=nrow(x)
          if (is.null(p1f)) {
            S0fns = function (x,t,p0f,ltm) {
              #this function calculates the survival curve for a fixed value of the covariate vector x
              p0=p0f(x,t)
              S0=rep(0,ltm)
              S0[1]=1-p0[1]
              for (i in 2:ltm) S0[i]=S0[i-1]*(1-p0[i])
              return(S0)
            }
            S0MSD=rowMSD(apply(x,1,S0fns,t=time,p0f=p0f,ltm=ltm))
            LTSlist=list(time=time,S0=S0MSD$rm,S0err=S0MSD$rsd/sqrt(n),accept.rate=acc.rate)
          } else {
            S01fns = function (x,t,p0f,p1f,ltm) {
              #this function calculates the survival curves for a fixed value of the covariate x
              p0=p0f(x,t)
              p1=p1f(x,t)
              S0=rep(0,ltm)
              S1=S0
              S0[1]=1-p0[1]
              S1[1]=1-p1[1]
              for (i in 2:ltm) {
                S0[i]=S0[i-1]*(1-p0[i])
                S1[i]=S1[i-1]*(1-p1[i])}
              return(c(S0,S1))
            }
            S01MSD=rowMSD(apply(x,1,S01fns,t=time,p0f=p0f,p1f=p1f,ltm=ltm))
            LTSlist=list(time=time,S0=S01MSD$rm[1:ltm],S0err=S01MSD$rsd[1:ltm]/sqrt(n),
                         S1=S01MSD$rm[(ltm+1):(2*ltm)],
                         S1err=S01MSD$rsd[(ltm+1):(2*ltm)]/sqrt(n),
                         accept.rate=acc.rate)
          }
        } else {
          sdV=sqrt(diag(V))
          #the integration is done after the covariates are transformed to lie in a subset of the unit hypercube by using the probit transformation
          L2=pnorm((L-M)/sdV)
          U2=pnorm((U-M)/sdV)
          if (is.null(p1f)) {
            S0fna = function (u,t,M,corV,sdV,p0f,ltm) {
              #this function calculates the survival curves for a fixed value of the coordinatewise probit-transformed covariate vector u.
              v=qnorm(u)
              p0=p0f(v*sdV+M,t)
              S0=rep(0,ltm)
              S0[1]=1-p0[1]
              for (i in 2:ltm) S0[i]=S0[i-1]*(1-p0[i])
              return(S0*dmvnorm(v,sigma=corV)*exp(sum(v^2)/2))
            }
            i1=adaptIntegrate(S0fna,L2,U2,t=time,M=M,corV=cov2cor(V),sdV=sdV,p0f=p0f,ltm=ltm,tol=1e-05,fDim=ltm,maxEval=100000,absError=0)
            S0=i1$integral*(2*pi)^(length(sdV)/2)/acc.rate
            S0err=i1$error*S0
            LTSlist=list(time=time,S0=S0,S0err=S0err,accept.rate=acc.rate)
          } else {
            S01fna = function (u,t,M,corV,sdV,p0f,p1f,ltm) {
              #this function calculates the survival curves for a fixed value of the coordinatewise probit-transformed covariate vector u.
              v=qnorm(u)
              p0=p0f(v*sdV+M,t)
              p1=p1f(v*sdV+M,t)
              S0=rep(0,ltm)
              S1=S0
              S0[1]=1-p0[1]
              S1[1]=1-p1[1]
              for (i in 2:ltm) {
                S0[i]=S0[i]*(1-p0[i])
                S1[i]=S1[i]*(1-p1[i])}
              return(c(S0,S1)*dmvnorm(v,sigma=corV)*exp(sum(v^2)/2))
            }
            i1=adaptIntegrate(S01fna,L2,U2,t=time,M=M,corV=cov2cor(V),sdV=sdV,p0f=p0f,p1f=p1f,ltm=ltm,tol=1e-05,fDim=2*ltm,maxEval=100000,absError=0)
            S01=i1$integral*(2*pi)^(length(sdV)/2)/acc.rate
            S01err=i1$error*S01
            LTSlist=list(time=time,S0=S01[1:ltm],S0err=S01err[1:ltm],
                         S1=S01[(ltm+1):(2*ltm)],S1err=S01err[(ltm+1):(2*ltm)],accept.rate=acc.rate)
          }
        }
      }
    return(LTSlist)
}

