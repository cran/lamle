\name{lamleout-class}
\Rdversion{1.1}
\docType{class}
\alias{lamleout-class}

\title{Class "lamleout"}
\description{
Output for an estimated generalized linear latent variable model.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("lamleout", ...)}.
}
\section{Slots}{
  \describe{
    \item{Data}{A list containing the supplied data (y, group) and information about the data (N, m).}
    \item{Estimates}{A list with parameter estimates (par, partable), standard errors (separ), mode estimates (map), covariance/information matrices (acov, Amat, Bmat), and various other quantities and objects derived from the parameter estimates (loadings, covmat, mu, modelpars, modelpartsoptC).}
    \item{Model}{A list containing model definitions (model, covstruct, modeltype), constraints (groupequal, parequal, parfix), estimation settings (method, accuracy) and filters.}
    \item{Optim}{A list containing the log-likelihood (loglik, logliki), gradient (grad, gradi), AIC, BIC, parameter estimate history (partrace), gradient history (glltrace), number of iterations (iter), and function call (call).}
    \item{Timing}{A data frame containing timing information of data preparation (prep), estimation (est) and post-estimation computations (se).}
  }
}
%%\section{Methods}{
%% \describe{
%%    \item{plot}{\code{signature(x = "lamleout")}: ... }
%%    \item{summary}{\code{signature(object = "lamleout")}: ... }
%%	 }
%%}

\author{
\email{bjorn.andersson@cemo.uio.no}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lamle}}
}
\examples{
showClass("lamleout")
}
\keyword{classes}
