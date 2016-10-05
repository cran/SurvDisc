\name{AsympDiscSurv}
\encoding{latin1}
\alias{AsympDiscSurv}
\title{Asymptotic Estimate of Mean and Variance of Log-hazard Ratio for Discrete Time Survival}
\description{
calculates the expected estimated log-hazard ratio and the estimated variance for large sample sizes when there are two groups with possibly non-proportional hazards and possible unequal randomization and censoring distributions.
}
\usage{
AsympDiscSurv(h0,h1,p0,p1,method=c("efron","breslow","PrenticeGloeckler"),tol=1E-12)}
\arguments{
  \item{h0}{vector of hazard rates in the control group}
  \item{h1}{vector of hazard rates in the treatment group}
  \item{p0}{vector of probabilities of being in the risk set and in the control group. See \code{Details} section below.}
  \item{p1}{vector of probabilities of being in the risk set and in the treatment group. See \code{Details} section below.}
  \item{method}{method for handling ties.}
  \item{tol}{a positive scalar giving the tolerance at which the maximum absolute value of the gradient is considered close enough to 0 to stop the algorithm.}
}

\details{
This calculates the asymptotic mean of the coefficient estimated by a proportional hazards regression model between two groups.

If there are \code{r} intervals, the vectors need only be of length r-1 since all subjects reaching the final interval will be assumed to have an event at some time in the last interval.

\code{p0} and \code{p1} are not the survival curves because they also include information about the allocation ratio between groups and the censoring distribution. The j^{th} element of \code{p0} is the probability of being assigned to the control group and being at risk at time \code{time[j]}. \code{p0+p1} is always less than or equal to 1 and should be close to 1 at the first time point and decreasing with time.  Note that subjects censored at \code{time[j]} are not in the risk set, only subjects who have an event at this time or later or who are censored later. This definition of censoring time is the definition used in the reference and may be different than used in other places. Add 1 to all censored times if desired to force censoring to conform with the more standard ways. With equal allocation and no censoring, then \code{p0[1]=p1[1]=0.5}.

}
\value{
  A list which contains:
   \item{coefficients}{the estimated coefficient (log-hazard ratio)}
   \item{varn}{the asymptotic variance multiplied by n where n is the total sample size combined in both groups}
   }

\author{John Lawrence,\email{john.lawrence@fda.hhs.gov}}

\seealso{\code{\link[SurvDisc]{LongToSurv}}}
\examples{
set.seed(1234)
nsim=1
n=250
k=50
trt=c(rep(0,n),rep(1,n))
betaef=rep(0,nsim)
varef=betaef
betapg=betaef
varpg=betaef
m1=3.05
for (i in 1:nsim){
  x=rexp(2*n,1)
  x[(n+1):(2*n)]=x[(n+1):(2*n)]/2
  x=ceiling(x*(k-1)/m1)
  x=pmin(x,k)
  cens=rbinom(2*n,1,0.9)
  pg1=PrenticeGloeckler.test(x,cens,trt,k)
  betapg[i]=pg1$coefficient
  varpg[i]=(pg1$coefficient/pg1$wald.test)^2
  efron=survival::coxph(survival::Surv(x,cens)~trt,ties="efron")
  betaef[i]=efron$coef
  varef[i]=efron$var}

h0=0.9*(exp(-c(0:(k-2))*m1/(k-1))-exp(-c(1:(k-1))*m1/(k-1)))
h0=h0/(h0+exp(-c(1:(k-1))*m1/(k-1)))
p0=exp(-c(0:(k-1))*m1/(k-1))
p0=(p0[1:(k-1)]*0.9+p0[2:k]*0.1)/2
h1=0.9*(exp(-c(0:(k-2))*2*m1/(k-1))-exp(-c(1:(k-1))*2*m1/(k-1)))
h1=h1/(h1+exp(-c(1:(k-1))*2*m1/(k-1)))
p1=exp(-2*c(0:(k-1))*m1/(k-1))
p1=(p1[1:(k-1)]*0.9+p1[2:k]*0.1)/2

fa=AsympDiscSurv(h0=h0, h1=h1,p0=p0,p1=p1)
c(fa$estimate,fa$var/(2*n))
c(mean(betaef),var(betaef),mean(varef))

}
