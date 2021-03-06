AsympDiscSurv <- function(h0,h1,p0,p1,method=c("efron","breslow","PrenticeGloeckler"),tol=1E-12)
{
  lenh0=length(h0)
  beta1=mean(h1/h0)
  beta0=beta1-1
  grad=1
  if (!is.element(substr(method[1],1,1),c("b","P"))) {
    while (abs(grad)>tol) {
      beta0=beta1
      theta=exp(beta0)
      grad=0
      hess=0
      for (i in 1:lenh0) {
        p0byp1th=p0[i]/(p1[i]*theta)
        est1=(h1[i]*(h1[i] + h0[i]*p0byp1th) + (h0[i] - h1[i])*p0byp1th*
                (log(1 + p0byp1th) - log(1 - h1[i] + p0byp1th - h0[i]*p0byp1th)))/(h1[i] + h0[i]*p0byp1th)^2

        est2=(-h1[i]^3 - h0[i]*h1[i]^2*p0byp1th +
                ((h0[i] - h1[i])^2*p0byp1th^2)/(-1 + h1[i] + (-1 + h0[i])*p0byp1th) +
                p0byp1th*(((h0[i] - h1[i])^2*p0byp1th)/(1 + p0byp1th) +
                            2*h1[i]*(-h0[i] + h1[i])*log(1 + p0byp1th)) +
                (h1[i] + h0[i]*p0byp1th)*(h1[i]*(h1[i] + h0[i]*p0byp1th) +
                                            (h0[i] - h1[i])*p0byp1th*(log(1 + p0byp1th) -
                                                                        log(1 - h1[i] + p0byp1th - h0[i]*p0byp1th))) -
                2*h1[i]*(-h0[i] + h1[i])*p0byp1th*log(1 - h1[i] + p0byp1th - h0[i]*p0byp1th))/(h1[i] + h0[i]*p0byp1th)^3
        grad=grad-h1[i]*p1[i]+est1*(h1[i]*p1[i]+h0[i]*p0[i])
        hess=hess+est2*(h1[i]*p1[i]+h0[i]*p0[i])}
      beta1=beta0-grad/hess}
    varn=1/hess} else if (substr(method[1],1,1)=="b") {
      d1=p1*h1
      d0=p0*h0
      while (abs(grad)>tol) {
        beta0=beta1
        theta=exp(beta0)
        grad=sum(-d1+(d0+d1)*(theta*p1)/(theta*p1+p0))
        hess=sum((d0+d1)*((theta*p1)/(theta*p1+p0)-(theta*p1)^2/(theta*p1+p0)^2))
        beta1=beta0-grad/hess}
      varn=1/hess} else {
        lpm1=length(p0)
        D0=p0*h0
        D1=p1*h1
        R0=p0
        R1=p1
        Ri=R0+R1
        N0i=D0+R0
        N1i=D1+R1
        Ni=N0i+N1i
        logNR=log(Ni/Ri)
        p1=c(log(logNR),0)
        p0=p1-1
        while (max(abs(grad))>tol) {
          p0=p1
          gam=p0[1:lpm1]
          beta=p0[lpm1+1]
          grad=exp(gam)*(D0/(1-exp(exp(gam)))+D1*exp(beta)/(1-exp(exp(gam+beta)))+R0+R1*exp(beta))
          grad=c(grad,sum(exp(beta+gam)*(D1/(1-exp(exp(gam+beta)))+R1)))
          hess=exp(gam)*(D0*(1+exp(exp(gam))*(exp(gam)-1))/(1-exp(exp(gam)))^2+
                           D1*exp(beta)*(1+exp(exp(beta+gam))*(exp(beta+gam)-1))/(1-exp(exp(beta+gam)))^2+R0+R1*exp(beta))
          hess=diag(c(hess,0))
          hess[lpm1+1,1:lpm1]=(D1*(1+exp(exp(beta+gam))*(exp(beta+gam)-1))+R1*(1-exp(exp(beta+gam)))^2)*exp(beta+gam)/(1-exp(exp(beta+gam)))^2
          hess[1:lpm1,lpm1+1]=hess[lpm1+1,1:lpm1]
          hess[lpm1+1,lpm1+1]=sum(hess[lpm1+1,1:lpm1])
          p1=p0-ginv(hess)%*%grad}
        beta0=p1[lpm1+1]
        varn=ginv(hess)[lpm1+1,lpm1+1]}
  return(list(estimate=beta0,varn=varn))}

