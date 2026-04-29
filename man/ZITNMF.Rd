\name{ZITNMF}
\docType{methods}
\alias{ZITNMF}

\title{
Zero-Inflated Tweedie Non-negative Matrix Factorization (ZITNMF)
}
\description{
The input data is assumed to be a non-negative matrix with excess zeros.
ZITNMF decomposes the matrix into two low-dimensional factor matrices
while accounting for zero inflation using a Tweedie distribution model.
}
\usage{
ZITNMF(X, Z=NULL, pseudocount=1e-10,
  initF=NULL, initA=NULL, initW=0.5, fixF=FALSE, fixA=FALSE,
  initializer=c("NMF", "ALS", "Random"),
  J=3, Beta=2, phi=1, thr=1e-10, num.iter=100, verbose=FALSE)
}

\arguments{
\item{X}{
The input non-negative matrix which has N-rows and M-columns.
}
\item{Z}{
The zero-inflation indicator matrix which has N-rows and M-columns
(Default: NULL). If NULL, it is estimated from the data.
}
\item{pseudocount}{
The pseudo count to avoid zero division (Default: 1e-10).
}
\item{initF}{
The initial values of factor matrix F, which has N-rows and J-columns
(Default: NULL).
}
\item{initA}{
The initial values of factor matrix A, which has M-rows and J-columns
(Default: NULL).
}
\item{initW}{
The initial value of zero-inflation probability w, a scalar between 0 and 1
(Default: 0.5).
}
\item{fixF}{
Whether the factor matrix F is updated in each iteration step
(Default: FALSE).
}
\item{fixA}{
Whether the factor matrix A is updated in each iteration step
(Default: FALSE).
}
\item{initializer}{
The initialization method for factor matrices
(Default: "NMF"). "NMF" uses Beta-divergence NMF, "ALS" uses SVD-based
initialization, and "Random" uses random uniform values.
}
\item{J}{
The number of low-dimension (Default: 3).
}
\item{Beta}{
The Beta parameter for Beta-divergence. Beta=2 corresponds to Frobenius norm,
Beta=1 to KL divergence, Beta=0 to IS divergence (Default: 2).
Note that Beta < 1 requires the \code{tweedie} package.
}
\item{phi}{
The dispersion parameter for the Tweedie distribution (Default: 1).
}
\item{thr}{
The threshold for convergence based on relative change of Beta-divergence
(Default: 1e-10). If NULL, no convergence check is performed and
all num.iter iterations are run.
}
\item{num.iter}{
The number of maximum iterations (Default: 100).
}
\item{verbose}{
If verbose is TRUE, the iteration number and Beta-divergence value
are printed (Default: FALSE).
}
}

\value{
  F : The factor matrix (N-rows and J-columns)

  A : The factor matrix (M-rows and J-columns)

  Z : The estimated zero-inflation indicator matrix (N-rows and M-columns)

  w : The estimated zero-inflation probability (scalar)

  RecError : The Beta-divergence values at each iteration

  RelChange : The relative change of Beta-divergence at each iteration
}

\references{
Zero-inflated models for count data and their extensions.
}

\author{Koki Tsuyuzaki}

\examples{
if(interactive()){
  # Test data
  X <- toyModel("NMF")

  # Frobenius
  out <- ZITNMF(X, J=3, Beta=2, num.iter=10, initializer="Random")

  # KL
  out <- ZITNMF(X, J=3, Beta=1, num.iter=10, initializer="Random")
}
}

\keyword{methods}
