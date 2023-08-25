\name{lamle.fit}
\alias{lamle.fit}
\title{
  Model Fit Statistics for an Estimated Latent Variable Model
}
\description{
A function to compute model fit statistics for latent variable models estimated with lamle. Currently only computes the standardized root mean square residuals (SRMSR) from the model-implied correlation matrix based on the estimated model parameters and the sample correlation matrix.
}
\usage{
lamle.fit(obj, obs, N = NULL, z = NULL, seed = NULL, 
          use = "complete", residmat = FALSE)

}
\arguments{
  \item{obj}{
An estimated model from function lamle.
}

  \item{obs}{
The observed data matrix, with missing values coded as NA. Ordinal data must be coded as successive integers starting from 1.
}

  \item{N}{
A numeric scalar indicating the multiplier of the sample size in the data to use when computing the model fit statistics.
}

\item{z}{
An optional matrix of latent variable values to use when computing the model fit statistics.
}

  \item{seed}{
The seed to use when simulating data.
}

  \item{use}{
A character vector indicating how the correlaion matrix used for the SRMSR statistic should be estimated. See function cor() for details.
}

  \item{residmat}{
Optional argument to return the full residual correlation matrix.
}

}
\value{
A numeric vector with fit statistics.
}

\author{
Björn Andersson <bjoern.h.andersson@gmail.com> and Shaobo Jin  <shaobo.jin@statistik.uu.se>
}

\references{
Andersson, B., and Xin, T. (2021). Estimation of latent regression item response theory models using a second-order Laplace approximation. \emph{Journal of Educational and Behavioral Statistics}, \bold{46(2)} 244-265. <doi:10.3102/107699862094519>

Jin, S., Noh, M., and Lee, Y. (2018). H-Likelihood Approach to Factor Analysis for Ordinal Data. \emph{Structural Equation Modeling: A Multidisciplinary Journal}, \bold{25(4)}, 530-540. <doi:10.1080/10705511.2017.1403287>

Shun, Z., and McCullagh, P. (1995). Laplace approximation of high dimensional integrals. \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, \bold{57(4)}, 749-760. <doi:10.1111/j.2517-6161.1995.tb02060.x>

}

\examples{
##### Load required package.
library(mvtnorm)

#####  Generate ordinal data from the GPCM
#####
##### Item parameter generation
set.seed(123)
GPCMa <- runif(10, 0.8, 2)
GPCMb <- vector("list", 10)
for(i in 1:10) GPCMb[[i]] <- -c(runif(1, -3, -2), 
                                runif(1, -1.5, -0.5), 
                                runif(1, 0, 1),
                                runif(1, 1.5, 2.5))
GMCMa2d <- matrix(0, nrow = 10, ncol = 2)
GMCMa2d[1:5, 1] <- GPCMa[1:5]
GMCMa2d[6:10, 2] <- GPCMa[6:10]

##### Latent variables in two groups
n <- 200
set.seed(1234)
covmat2d <- diag(rep(1, 2))
covmat2d[!diag(2)] <- 0.6
latmat2dA <- rmvnorm(n, c(0, 0), covmat2d)
latmat2dB <- rmvnorm(n, c(1, 1), covmat2d)

##### Observed data
set.seed(12345)
dataGPCM2d <- matrix(NA, nrow = 2 * n, ncol = 10)
dataGPCM2d[1:n, ] <- DGP(GMCMa2d, GPCMb, rep("GPCM", 10), 
                         latmat2dA)
dataGPCM2d[(n + 1):(2 * n), ] <- DGP(GMCMa2d, GPCMb, 
                                     rep("GPCM", 10), latmat2dB)

#####  Setup two-dimensional independent-clusters model
mydim2d <- matrix(NA, nrow = 10, ncol = 2)
mydim2d[1:5, 1] <- 1
mydim2d[6:10, 2] <- 1

#####  Lap(2)
estLap2 <- lamle(y = dataGPCM2d[1:n,], model = mydim2d, 
                 modeltype = rep(c("GPCM", "GRM"), 5), 
                 method = "lap", accuracy = 2, optimizer = 
                 "BFGS", inithess = "crossprod")

lamle.fit(estLap2, dataGPCM2d[1:n,], N = 10)
}