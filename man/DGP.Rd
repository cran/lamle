\name{DGP}
\alias{DGP}
\title{
  Generation of Observed Data From a Generalized Linear Latent Variable Model
}
\description{
A function to generate simulated observed data based on a generalized linear latent variable model, supporting ordinal, count, and continuous observed variables.
}
\usage{
DGP(a, b, modeltype, z)
}
\arguments{
  \item{a}{
  A matrix with number of rows equal to the number of observed variables and number of columns equal to the number of latent variables, containing the discrimination parameters for each item and dimension.
}

  \item{b}{
  A list with length equal to the number of observed variables, containing intercept parameters and, if applicable, the scale parameter for each observed variable.
}

  \item{modeltype}{
  A character vector with length equal to the number of observed variables to simulate, indicating which model should be used for each variable. Options are "GPCM", "GRM", "negbin", "normal", and "poisson".
}

  \item{z}{
  The matrix of latent variable values, with rows indicating the individual and the columns indicating the dimension.
}

}
\value{
A data matrix with the observed item responses.
}

\author{
Björn Andersson <bjoern.h.andersson@gmail.com>
}
\references{
Muraki, E. (1992). A Generalized Partial Credit Model: Application of an EM Algorithm. \emph{Applied Psychological Measurement}, \bold{16}, 159-176.
}

\examples{
##### Load required package.
library(mvtnorm)

##### Item parameters
set.seed(123)
GPCMa <- runif(60, 0.8, 2)
GPCMb <- vector("list", 60)
for(i in 1:60) GPCMb[[i]] <- -c(runif(1, -3, -2), 
                                runif(1, -1.5, -0.5), 
                                runif(1, 0, 1),
                                runif(1, 1.5, 2.5))
								
##### Generate data from a two-dimensional independent-clusters model
mydata <- DGP(matrix(c(GPCMa[1:30], rep(0, 60), GPCMa[31:60]), 
              ncol = 2), GPCMb, model = rep("GPCM", 60), 
              rmvnorm(1000, c(0, 0), matrix(c(1, 0.6, 0.6, 1), 
              ncol = 2)))
}