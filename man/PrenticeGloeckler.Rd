\name{PrenticeGloeckler.test}
\alias{PrenticeGloeckler.test}
\title{Regression for Grouped Survival Data Function}
\description{
This function calculates the estimated hazard ratio for grouped survival data described in the reference below.
}
\usage{
PrenticeGloeckler.test(time,event,grp,r)}
\arguments{
  \item{time}{vector of times to event or censoring. The times are assumed to be integers from 1, 2, .., r corresponding to the discrete time points or the continuous time intervals A1, ..., Ar}
  \item{event}{vector of binary status indicator variables (0 = censored at the start of the interval, 1 = event during the interval)}
  \item{grp}{vector of binary group indicators (0 or 1)}
  \item{r}{number of time points or intervals}
}

\details{
The hazard functions and hazard ratio are estimated for grouped survival data.
}
\value{
	A list consisting of:
    \item{coefficient}{The estimated coefficient (log hazard ratio) found by maximizing the likelihood.}
    \item{indx}{vector of time points where the hazard functions are estimated. The subset of \code{1,...,r-1} with at least one event.}
    \item{gamma}{numeric vector with the same length as \code{indx} representing the log(-log(hazard rate)) in the control group for time points in the vector \code{indx}}
    \item{grad1}{gradient evaluated at \code{(gamma[indx],ceofficient)}}
    \item{r}{number of time points or time intervals}
    \item{hess1}{hessian matrix evaluated at the maximum likelihood estimate.}
    \item{ll0}{log-likelihood evaluated at ceofficient=0. includes attributes \code{"gradient"} and \code{"hessian"}}
    \item{ll1}{log-likelihood at maximum likeohood estimate.  includes attributes \code{"gradient"} and \code{"hessian"}}
    \item{score.test}{value of the score test statistic for testing coefficient=0 (see reference).}
    \item{lr.test}{value of the likelihood ratio test statistic, 2*(ll0-ll1)}
    \item{wald.test}{value of the Wald test statistic; the estimated coefficient divided by the square root of the estimated variance.}
}

\references{
Prentice, R. L. and Gloeckler, L.A. (1978). Regression analysis of grouped survival data with application to breast cancer data. \emph{Biometrics}, 57 -- 67
}

\author{John Lawrence}
\examples{
set.seed(1234)
nsim=1
n=250
tn=2*n
k=0.1*tn
betaef=rep(0,nsim)
betapg=rep(0,nsim)
cens=rep(1,2*n)
trt=c(rep(0,n),rep(1,n))

for (i in 1:nsim) {
  x=rexp(tn,1)
  x[(n+1):tn]=x[(n+1):tn]/2
  m1=max(x[(n+1):tn])
  x=ceiling(x*(k-1)/m1)
  x[(n+1):tn]=pmin(x[(n+1):tn],k-1)
  x[1:n]=pmin(x[1:n],k)
  pg1=PrenticeGloeckler.test(x,cens,trt,k)
  betapg[i]=pg1$coefficient
  betaef[i]=survival::coxph(survival::Surv(x,cens)~trt,ties="efron")$coef}
mean(betaef)
mean(betapg)
}
