\name{simGFRdata}
\docType{data}
\alias{simGFRdata}
\alias{simGFR}
\alias{fix.beta}
\alias{fixed.leGFR}
\alias{fixed.leGFRtrt}
\alias{fixed.time}
\alias{fixed.trt}
\alias{fixed.trttime}
\alias{mu.base.leGFR}
\alias{var.base.leGFR}
\alias{res.sd}
\alias{sig.intercept}
\alias{sig.time}
\alias{sig.cor}
\title{Data Set Containing Simulated Longitudinal eGFR}
\description{A data set containing simulated values of log-eGFR measured longitudinally over time as a function of baseline eGFR.
The data were simulated from a mixed effects model with the following form (using the lme model structure syntax; see format section below
for definition of variables):

cfb ~ time + x1:time + trt + trt:time + trt:x1:time + 0

and these coefficients:

       time         trt     time:x1    time:trt time:x1:trt
-0.6447911  -0.0478315   0.1333391   0.2186963  -0.0458998

In addition, each subject has a random slope and intercept. The baseline eGFR were simulated from a log-Normal distribution.
}

\usage{data(simGFRdata)}

\format{
  Fixed effect coefficients used to simulate the data:
  fix.beta
  fixed.leGFR
  fixed.leGFRtrt
  fixed.time
  fixed.trt
  fixed.trttime

  mu.base.leGFR: mean of baseline log-eGFR
  var.base.leGFR: variance of baseline log-eGFR

  res.sd: residual error standard deviation. note this is for a single log-eGFR, so the standard deviation for the change from baseline is sqrt(2)*res.sd and the residual error for cfb within a patient have correlation 0.5.

  Variance components of random effects distribution:
  sig.intercept: standard deviaiton of random intercept
  sig.time: standrd deviation of random slope
  sig.cor: correlation

	A data frame named simGFR that consists of fourteen columns and 28800 rows. The variables are:
	PID: patient ID
	trt: the treatment group indicator
	x1: measured value of baseline log-eGFR
	time: time from baseline measured in years
  alphai: subject's random intercept
	betai: subject's random slope
	alpha: subject's intercept including fixed and random effects
	beta: subject's slope including fixed and random effects
	cfb0: the measurement error for the baseline log-eGFR
  x: the unobserved "true" baseline log-eGFR
  cfb: the change from baseline in measured log-eGFR
}

\keyword{datasets}
