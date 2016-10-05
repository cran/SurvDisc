\name{TKVsurv}
\docType{data}
\alias{TKVsurv}
\title{Data Set Containing Fitted Model Estimates and Covariances}
\description{A data set containing the estimates fromthe fitted Cox proportional hazards model from a dataset of patients with Autosomal Dominant Polycystic Kideny Disease. See references for further details. The model has 6 parameters describing how the hazard changes for different levels of the 3 covariates. In addition, there are 3 strata correspoding to the different imaging modalities: CT, MRI, US.
}

\usage{data(TKVsurv)}

\format{
  A list with the following components
	covariates: names of covariates
	mean: sample mean of covariates
	stand.dev.: standard deviation of covariates
	labels: labels for coefficients and rows and columns of covariance matrix sigma
	beta: estimated coefficients in the proportional hazards model
	sigma: estimated covariance martix for beta
	CT.time: time points where Survival curve changes for CT strata
	CT.lcumhaz: estimated log-cumulative hazard in CT strata with coefficients = 0
	CT.sig17: estimated 7 elements to fill the last row (and last column) of the covariance matrix
  MRI.time: time points where Survival curve changes for MRI strata
	MRI.lcumhaz: estimated log-cumulative hazard in MRI strata with coefficients = 0
	MRI.sig17: estimated 7 elements to fill the last row (and last column) of the covariance matrix
	US.time: time points where Survival curve changes for US strata
	US.lcumhaz: estimated log-cumulative hazard in US strata with coefficients = 0
	US.sig17: estimated 7 elements to fill the last row (and last column) of the covariance matrix
}

\references{
http://www.fda.gov/downloads/Drugs/DevelopmentApprovalProcess/DrugDevelopmentToolsQualificationProgram/UCM458523.pdf
}
\author{John Lawrence, Jianjin Xu, Jim Hung, Sue Jane Wang}
\keyword{datasets}
