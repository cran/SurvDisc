\name{simexlme}
\encoding{latin1}
\alias{simexlme}

\title{SIMEX algorithm for linear mixed effects models}
\description{
Implementation of the SIMEX algorithm for measurement error models according
  to Cook and Stefanski.
}
\usage{
simexlme(model, model.model, SIMEXvariable, respvar, grpvar, corform, measurement.error,
  measurement.error.resp, lambda = c(0.5, 1, 1.5, 2), B = 100,
  fitting.method = "quadratic", jackknife.estimation = "quadratic")
}
\arguments{
  \item{model}{naive model}
  \item{model.model}{dataframe containing all variables in the model}
  \item{SIMEXvariable}{character name of the variable with measurement error. Assumed to be the baseline measurement.}
  \item{respvar}{character name of the response variable. The response is assumed to represent a change from baseline.}
  \item{grpvar}{character name of the grouping variable for the random effects in the model.}
  \item{corform}{formula for the correlation of residual errors within groups. see example}
  \item{measurement.error}{The known standard deviation of measurement errors for \code{SIMEXvariable}.}
  \item{measurement.error.resp}{The known stadard deviaiton for \code{respvar}}
  \item{lambda}{vector of lambdas for which the simulation step should be done}
  \item{B}{number of iterations for each lambda}
  \item{fitting.method}{fitting method for extrapolation. Only \code{linear} or \code{quadratic} are recommended.}
  \item{jackknife.estimation}{specifying the extrapolation method for jackknife variance estimation.}
}

\details{
See documentation for \code{\link[simex]{mcsimex}} function. This function for lme models was adapted from that function, which is designed to handle linear and generalized linear models, but not lme models. In this function, the measurement error variable must be the baseline value of some measurement and the response is the change from baseline in the same measurement. There is assumed to be one value of this baseline measurement per level of the grouping variable in the mixed effect model. The correlation between the measurement errors for two response values within a subject is assumed to be equal to be equal to the variance of baseline divided by the sum of the variance of baseline and variance of post-baseline errors. For example, for a study measuring the effect of some weight loss treatment, the grouping variable could be subject, the baseline weight is the covariate with measurement error and the response is change from baseline in weight.}
\value{
  An object of class 'simex' which contains:
   \item{coefficients}{the corrected coefficients of the SIMEX model}
   \item{SIMEX.estimates}{the estimates for every lambda}
   \item{model}{the naive model}
   \item{measurement.error}{the known error standard deviations for \code{SIMEXvariable}}
   \item{B}{the number of iterations}
   \item{extrapolation}{the model object of the extrapolation step}
   \item{fitting.method}{the fitting method used in the extrapolation step}
   \item{residuals}{the residuals of the main model}
   \item{fitted.values}{the fitted values of the main model}
   \item{call}{the function call}
   \item{variance.jackknife}{the jackknife variance estimate}
   \item{extrapolation.variance}{the model object of the variance extrapolation}
   \item{variance.jackknife.lambda}{the data set for the extrapolation}
   \item{variance.asymptotic}{the asymptotic variance estimates}
   \item{theta}{the estimates for every B and lambda}
   }

\references{
Cook, J.R. and Stefanski, L.A. (1994) Simulation-extrapolation estimation in
  parametric measurement error models. \emph{Journal of American Statistical
  Association}, \bold{89}, 1314 -- 1328
}

\author{John Lawrence,\email{john.lawrence@fda.hhs.gov}, Jianjin Xu, Wolfgang Lederer, Heidi Seibold}

\seealso{\code{\link[simex]{simex}},\code{\link[nlme]{lme}}}
\examples{
set.seed(1234)
data("simGFRdata")
simGFR=simGFR[is.element(simGFR$time,c(1:12)/4) & is.element(simGFR$PID,c(1:80)*100),]

fm2=nlme::lme.formula(fixed = cfb ~ time + x1:time + trt + trt:time + trt:x1:time + 0,
                data = simGFR, random = ~time | PID,
                correlation = nlme::corCompSymm(0.5,form = ~time | PID, fixed = TRUE),
                control=nlme::lmeControl(returnObject=TRUE))

(s1 = simexlme(model=fm2, model.model=simGFR[,c("cfb","PID","time","x1","trt")],
                SIMEXvariable="x1",respvar="cfb",grpvar="PID",corform="~time | PID",
                measurement.error=res.sd,measurement.error.resp=res.sd,
                lambda = c(0.5,2),B = 2, fitting.method = "linear",
                jackknife.estimation = FALSE))

plot(s1)

#values of fixed effects used to simulate data
c(fixed.time,fixed.trt,fixed.leGFR,fixed.trttime,fixed.leGFRtrt)

#naive estimates
fm2$coefficients$fixed

#SIMEX corrected estimates
s1$coefficients
}
