\name{rowMSD}
\alias{rowMSD}
\title{Mean and Standard Deviation estimates for each row in a matrix}
\description{
  Calculates the sample mean and standard deviation for each row in a matrix. The mean vector is calculated first. The elements of the matrix are then centered by the mean vector before the sample standard deviation is calculated.
}
\usage{
rowMSD(x)
}
\arguments{
  \item{x}{numeric matrix}
}
\value{
	A list consisting of:
    \item{rm}{vector of row means}
    \item{rsd}{vector of row standard deviations}
}
\author{John Lawrence}
\examples{
x=matrix(rnorm(1000),nrow=10)
rowMSD(x)
}


