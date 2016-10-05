\name{LongToSurv}
\alias{LongToSurv}
\title{Longitudinal To Survival Function}
\description{
This function calculates the survival curve for events where the events are defined by some function of a variable measured longitudinally. The events can be defined with or without confirmation (see arguments and details).  The survival curve is integrated over a distribution of covariates. The term "covariates" is used loosely here and includes all terms in the mixed effects longitudinal model including random effects and error terms. This distribution is assumed to be truncated multivariate normal.
}
\usage{
LongToSurv(M,V,L,U,time,p0f,p1f=NULL,method=c("simulation","asymptotic"),
  conf.type=c("scheduled","unscheduled","none"),nsim=100000)}
\arguments{
  \item{M}{Mean vector for the parent multivariate normal distribution of the covariates.}
  \item{V}{Covariance matrix for the parent multivariate normal distribution of the covariates.}
  \item{L}{vector of lower limits for the covariates.}
  \item{U}{vector of upper limits for the covariates.}
  \item{time}{vector of time points.}
  \item{p0f}{multi-valued function that calculates the probability of crossing the threshold at each sceduled visit time point in the control group. If \code{method="unscheduled"}, the probability of crossing the threshold at both the scheduled and subsequent unscheduled visit; even the last visit is assumed to allow an unscheduled confirmation visit.}
  \item{p1f}{Optional multi-valued function that calculates the probability of crossing the threshold in the treatment group.}
  \item{method}{Method used; either \code{"simulation"} or \code{"analytic"}. Defaults to \code{"simulation"}.}
  \item{conf.type}{type of confirmation. "none" meaing a single value crossing the threshold is an event, "scheduled" meaning two consecutive scheduled measurements crossing the threshold, or "unscheduled" meaning that after a qualifying event at a scheduled visit, a subsequent measurement is taken at an unscheduled visit to potentially confirm the event.}
  \item{nsim}{Approximate number of simulated covariate values used. Used only if method = \code{"simulation"}.}
}

\details{
The discrete survival function is found given a distribution of covariates and a longitudinal model. The event is defined by the response variable crossing a threshold value either once (confirmation = "none") or twice in successive time points.  The distribution of the covariates is assumed to be truncated multivariate normal. If method is \code{"simulation"}, then /code{nsim/accept.rate} values of the covariates are simulated first. The truncation conditions are tested and approximately \code{nsim} of these covariates will be accepted. The survival curve is found and averaged over the covariate values in the sample. If the method is \code{"analytic"}, then the survival curve function is integrated analytically (using the \code{adaptIntegrate} function from the \code{cubature} package).
}
\value{
	A list consisting of:
    \item{times}{numeric vector of time points}
    \item{S0}{numeric vector of survival beyond time t in the control group}
    \item{S0err}{numeric vector of the estimated standard error (or estimated absolute error for analytic method) of S0}
    \item{S1}{numeric vector of survival beyond time t in the test group.}
    \item{S1err}{numeric vector of the estimated standard error (or estimated absolute error for analytic method) of S1}
    \item{accept.rate}{estimate probability that a covariate vector from the parent multivariate normal disitrbution will lie between the truncation limits L and U.}
}
\author{John Lawrence}
\examples{
mu.AGE = 38.582
mu.lbtkv = 6.9276
mu.base.leGFR = 4.2237
var.AGE = 220.73
var.lbtkv = 0.46848
var.base.leGFR=0.19770
cov.AGE.lbtkv = 3.4075
cov.AGE.leGFR = -4.5065
cov.lbtkv.leGFR = -0.16303
sig.intercept=0.03975
sig.time=0.04505
sig.cor=0.008
res.sd=0.11470307/sqrt(2)

M=c(mu.AGE,mu.lbtkv,mu.base.leGFR,0,0,0)
V=diag(c(var.AGE,var.lbtkv,var.base.leGFR+res.sd^2,res.sd^2,sig.intercept^2,sig.time^2))
V[1,2] = V[2,1] = cov.AGE.lbtkv
V[1,3] = V[3,1] = cov.AGE.leGFR
V[2,3] = V[3,2] = cov.lbtkv.leGFR
V[3,4] = V[4,3] = V[4,4]
V[5,6] = V[6,5] = sig.cor*sig.intercept*sig.time
L=c(18,6.9,3.9,-Inf,-Inf,-Inf)
U=c(40,8,5,Inf,Inf,Inf)
time=c(1:12)/4

p0f=function(x,t) {
  fixed.time=-0.337166
  fixed.age=0.0008176
  fixed.lbtkv=-0.02409
  fixed.leGFR0=0.09591
  trt.acute=-0.047759
  trt.chronic=0.0191574
  res.sd=0.11470307/sqrt(2)
  pnorm((log(0.7)-as.vector(x[5]+outer(x[6]+fixed.age*x[1]+fixed.lbtkv*x[2]+
  fixed.leGFR0*(x[3]-x[4])+fixed.time,t)-x[4]))/res.sd)
}

p1f=function(x,t) {
  fixed.time=-0.337166
  fixed.age=0.0008176
  fixed.lbtkv=-0.02409
  fixed.leGFR0=0.09591
  trt.acute=-0.047759
  trt.chronic=0.0191574
  res.sd=0.11470307/sqrt(2)
  pnorm((log(0.7)-as.vector(x[5]+trt.acute+outer(x[6]+fixed.age*x[1]+fixed.lbtkv*x[2]+
  fixed.leGFR0*(x[3]-x[4])+fixed.time+trt.chronic,t)-x[4]))/res.sd)
}

LTS1=LongToSurv(M,V,L,U,time,p0f,p1f,nsim=100) #nsim much larger than 100 is recommended
LTS1
#LTS2=LongToSurv(M,V,L,U,time,p0f,p1f,method="analytic")
#LTS2
}
