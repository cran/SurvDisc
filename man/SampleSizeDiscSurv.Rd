\name{SampleSizeDiscSurv}
\alias{SampleSizeDiscSurv}
\title{Sample Size for Discrete Time Survival}
\description{
Calculates the sample size needed to achieve any given power for any specifed type 1 error rate.
}
\usage{
SampleSizeDiscSurv(power=0.9,alpha=0.025,alternative=c("less","greater"),beta0=0,
  h0,h1,p0,p1,ties.method=c("efron","breslow","PrenticeGloeckler"),
  method=c("asymptotic","simulation"),tol,AMV=NULL,nsim=10000,Nvec=NULL,
  test=c("Wald","Score"))
  }
\arguments{
  \item{power}{scalar value of the desired power. Default value is 0.9.}
  \item{alpha}{scalar value of the one-sided type 1 error rate. Default value is 0.025.}
  \item{alternative}{character specifying the type of alternative}
  \item{beta0}{scalar value of the log-hazard ratio on the boundary of the null hypothesis. Default is 0.}
  \item{h0}{vector of length r-1 containing the postulated hazard rates in the control group for the times 1, ..., r-1 or corresponding time intervals. Assumed to be r intervals with the last interval being infinite.}
  \item{h1}{vector of postulated hazard rates in the treatment group}
  \item{p0}{vector of probabilities of being in the risk set and in the control group. See \code{Details} section below.}
  \item{p1}{vector of probabilities of being in the risk set and in the treatment group. See \code{Details} section below.}
  \item{ties.method}{method for handling ties.}
  \item{method}{character specifiying the asymptotic or simalution based method for determining the sample size.}
  \item{tol}{a positive scalar giving the tolerance at which the maximum absolute value of the gradient is considered close enough to 0 to stop the algorithm used if \code{method="asymptotic"}.}
  \item{AMV}{AsympDiscSurv object from a previous call to the AsympDiscSurv function.}
  \item{nsim}{number of simulations used per \code{N} value in the \code{Nvec} vector. Used only if \code{method="simulation"}}
  \item{Nvec}{vector of sample sizes used in simulation based method. If none specifed, default is to use two N values close to the estimate from the asymptotic method (see details below).}
  \item{test}{character specifying the type of test statistics used. Used only for simulation based method because asymptotically, the tests are equivalent.}
}

\details{
If \code{method="asymptotic"}, then the mean of the test statistic (wald or score, which are equivalent asymptotically) for a sample size divided by sqrt(N) converges to a constant. This constant is found from the parameters in the result of the call to AsympDiscSurv. If the AsympDiscSurv object has already been found, it can be passed to this function in the arguments. If not, then this function calls AsympDiscSurv to find those paraemters.

If \code{method="simulation"}, then the mean of the test statistic is found for each sample size in the \code{Nvec} vector. The mean and variance of the test statistic for each \code{N} is found. Then, a linear regression is used to find the sample size that will provide the correct power. Each test statistic is asumed to have a mean that depends on sqrt(N) and the same variance. Theoretically, the variance should be close to 1, but the variance is estimated from the simulated values (not assumed equal to 1). The normality assumption is usually satisfied if the number of events is sufficiently large.

Neither the simulation nor the asymptotic method is reliable if the expected number of events is small (say, less than 20). The asymptotic method is faster. However, the simulation method has several advantages. First, the asymptotic variance found by the AsympDiscSurv function can differ from the true variance by a few percent even for moderately large sample sizes. The simulation based method estimates the true variance by simulation. Second, for moderatley large sample sizes, the score test can be different from the Wald test. Third, asymptotically the mean of the test statistic is approximately constant times sqrt(N), i.e. a linear function of sqrt(N) with no intercept. But, for small N, the relationship may not be so simple. The simulation method models the relationship for values of N close to the target value without making this strong assumption. The simulation method still assumes that the test statistic is normally distributed, so may be inaccurate for very small sample sizes or rare events.

Iy is assumed there are \code{r} time intervals, the vectors defining the hazard and at-risk rates have length \code{r-1} since all subjects reaching the final interval must have an event at some time in the last interval.

\code{p0} and \code{p1} are not the survival curves because they also include information about the allocation ratio between groups and the censoring distribution. The j^{th} element of \code{p0} is the probability of being assigned to the control group and being at risk at time \code{time[j]}. \code{p0+p1} is always less than or equal to 1 and should be close to 1 at the first time point and decreasing with time.  Note that subjects censored at \code{time[j]} are not in the risk set, only subjects who have an event at this time or later or who are censored later. This definition of censoring time is the definition used in the reference and may be different than used in other places. Add 1 to all censored times if desired to force censoring to conform with the more standard ways. With equal allocation and no censoring, then \code{p0[1]=p1[1]=0.5}.

}

\value{
  An object of class SSDS, which is a list containing:
  \item{N}{sample size that should provide the correct power}
  \item{alternative}{character specifying the type of alternative}
  \item{beta0}{scalar value of the log-hazard ratio on the boundary of the null hypothesis. Default is 0.}
  \item{ties.method}{method for handling ties.}
  \item{method}{character specifiying the asymptotic or simalution based method for determining the sample size.}
  \item{AMV}{AsympDiscSurv object}
  \item{EZobj}{required expected value of the test statistic}
  \item{Nvec}{vector of sample sizes used in the simulation}
  \item{EZvec}{vector of mean values of the test statistic for each value of N}
  \item{VZvec}{vector of sample variances for each value of N}
  \item{int.est}{estimate of the intercept in the linear relationship between sqrt(N) and expected value of the test statistic.}
  \item{slope.est}{estimate of the slope in the linear relationship between sqrt(N) and expected value of the test statistic.}
  \item{nsim}{number of simulations used per \code{N} value in the \code{Nvec} vector. Used only if \code{method="simulation"}}
  \item{test}{character specifying the type of test statistics used. Used only for simulation based method because asymptotically,
  the tests are equivalent.}
  }

\author{John Lawrence,\email{john.lawrence@fda.hhs.gov}}

\seealso{\code{\link[SurvDisc]{LongToSurv}},\code{\link[SurvDisc]{AsympDiscSurv}}}
\examples{
set.seed(1234)
k=50
m1=3.05

h0=0.9*(exp(-c(0:(k-2))*m1/(k-1))-exp(-c(1:(k-1))*m1/(k-1)))
h0=h0/(h0+exp(-c(1:(k-1))*m1/(k-1)))
p0=exp(-c(0:(k-1))*m1/(k-1))
p0=(p0[1:(k-1)]*0.9+p0[2:k]*0.1)/2
h1=0.9*(exp(-c(0:(k-2))*2*m1/(k-1))-exp(-c(1:(k-1))*2*m1/(k-1)))
h1=h1/(h1+exp(-c(1:(k-1))*2*m1/(k-1)))
p1=exp(-2*c(0:(k-1))*m1/(k-1))
p1=(p1[1:(k-1)]*0.9+p1[2:k]*0.1)/2

fa=AsympDiscSurv(h0=h0,h1=h1,p0=p0,p1=p1)

(SSDS1=SampleSizeDiscSurv(power=0.9,alpha=0.025,alternative="greater",beta0=0,h0,h1,
  p0,p1,ties.method="efron",method="asymptotic",AMV=fa,Nvec=NULL,test="Wald"))

}
