PrenticeGloeckler.test=function(time,event,grp,r) {
  #calculates the log-likelihood for grouped survival data with two groups from eq. (9) Prentice Gloeckler
  #time is a vector with the event times
  #the times are assumed to be integers from 1, 2, .., r, i.e. indicating the time interval A1, ..., Ar
  #event is a vector with the censoring indicators: 1 is an observed event, 0 is censored (at the beginning of the interval)
  #grp is a vector with the treatment group assignment: 0 for control and 1 for test group

  n0=sum(grp==0)
  n1=sum(grp==1)
  n=n0+n1
  if (max(abs(time - round(time)))>.Machine$double.eps^0.5 | min(time)<1 | max(time)>r) print("Times in time vector are not all positive integers less than or equal to r.") else {

    lbls=labels(table(time[event==1 & time<max(time)]))[[1]]
    l0=length(lbls)
    keep=c(1:n)[is.element(time,as.numeric(lbls))]
    grp=grp[keep]
    event=event[keep]
    time=time[keep]
    x0=time[grp==0]
    x1=time[grp==1]

    D0=table(x0)[lbls]
    D0[is.na(D0)]=0
    D1=table(x1)[lbls]
    D1[is.na(D1)]=0
    R0=n0-cumsum(D0)
    R1=n1-cumsum(D1)
    D0=table(time[grp==0 & event==1])[lbls]
    D0[is.na(D0)]=0
    D1=table(time[grp==1 & event==1])[lbls]
    D1[is.na(D1)]=0

    #D0[i] and D1[i] are the number of events at time i
    #R0[i] and R1[i] are the number of subjects surviving past time i

    ll=function(p,D0,D1,R0,R1,lpm1) {
      gam=p[1:lpm1]
      beta=p[lpm1+1]
      res=-sum(D0*log(1-exp(-exp(gam)))+D1*log(1-exp(-exp(gam+beta)))-
                 R0*exp(gam)-R1*exp(gam+beta))

      grad=exp(gam)*(D0/(1-exp(exp(gam)))+D1*exp(beta)/(1-exp(exp(gam+beta)))+R0+R1*exp(beta))
      grad=c(grad,sum(exp(beta+gam)*(D1/(1-exp(exp(gam+beta)))+R1)))
      attr(res,"gradient") =grad

      hess=exp(gam)*(D0*(1+exp(exp(gam))*(exp(gam)-1))/(1-exp(exp(gam)))^2+
                       D1*exp(beta)*(1+exp(exp(beta+gam))*(exp(beta+gam)-1))/(1-exp(exp(beta+gam)))^2+R0+R1*exp(beta))
      hess=diag(c(hess,0))
      hess[lpm1+1,1:lpm1]=(D1*(1+exp(exp(beta+gam))*(exp(beta+gam)-1))+R1*(1-exp(exp(beta+gam)))^2)*exp(beta+gam)/(1-exp(exp(beta+gam)))^2
      hess[1:lpm1,lpm1+1]=hess[lpm1+1,1:lpm1]
      hess[lpm1+1,lpm1+1]=sum(hess[lpm1+1,1:lpm1])
      attr(res,"hessian") =hess

      res}

    N0i=D0+R0
    N1i=D1+R1
    Ri=R0+R1
    Ni=N0i+N1i
    Di=D0+D1
    logNR=log(Ni/Ri)
    p=c(log(logNR),0)

    qi=(Ni*Ri/Di)*logNR^2
    score.test=sum(((D0*N1i-D1*N0i)/Di)*logNR)^2/sum(qi*N0i*N1i/Ni^2)

    ll0=ll(p,D0,D1,R0,R1,l0)

    p=nlm(ll,-ginv(attr(ll0,"hessian"))%*%attr(ll0,"gradient"),D0=D0,D1=D1,R0=R0,R1=R1,lpm1=l0,hessian=T)
    ll1=p$minimum
    attr(ll1,"gradient")=p$gradient
    attr(ll1,"hessian")=p$hessian

    l1=list(coefficient=p$est[l0+1],indx=lbls,gamma=p$est[1:l0],grad1=p$gradient,r=r,
            hess1=p$hessian,ll0=ll0,ll1=ll1,
            score.test=score.test,lr.test=2*as.numeric(ll0-ll1),
            wald.test=p$est[l0+1]/sqrt(ginv(p$hessian)[l0+1,l0+1]))
    class(l1)="PrenticeGloeckler.test"
    return(l1)
  }
}


print.PrenticeGloeckler.test=function(x, ...) print(x$coefficient)

