\name{lamle.predict}
\alias{lamle.predict}
\title{
  Compute Latent Variable Estimates from an Estimated Latent Variable Model
}
\description{
A function to plot different quantities for latent variable models estimated with lamle. Currently supports plotting of category probability functions, expected score functions, observed variable information functions and test information functions.
}
\usage{
lamle.predict(x, obs, estimator = "MAP", information = NULL, z.tol = 1e-8, group = NULL)

}
\arguments{
  \item{x}{
An estimated model from function lamle.
}

  \item{obs}{
A matrix of observations, with rows indicating cases and columns indicating variables. The number of columns and order of the columns must match the data used to estimate the model provided. If only partial data are available, specify NA for the missing data.
}

  \item{estimator}{
The estimator used for the scoring, with options "MAP" for the maximum aposteriori estimator (also known as Bayesian modal estimation or the posterior mode) and "MLE" for the maximum likelihood estimator. Note that the MLE is not defined if the responses are all in the lowest or highest categories for ordinal data. Default is "MAP".
}

  \item{information}{
Option to compute the information matrices evaluated at the latent variable estimates. Options are 'expected' for the expected information (Fisher information) and 'observed' for the observed information. Default is NULL, which avoids computing the information matrices.
}

  \item{z.tol}{
The tolerance used in the numerical optimization algorithm. Default is 1e-8.
}

  \item{group}{
A vector with length equal to the number of rows in argument 'obs' which indicates group membership of each case. Default is NULL, which sets all group memberships to the first group.
}

}
\value{
A list with latent variable estimates (a matrix) in the first entry and the corresponding information matrices (a list of matrices) in the second entry. 
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

lamle.predict(estLap2, dataGPCM2d[1:n,], estimator = "MAP")
}