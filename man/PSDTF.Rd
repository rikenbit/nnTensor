\name{PSDTF}
\docType{methods}
\alias{PSDTF}

\title{
Positive Semidefinite Tensor Factorization with LogDet Divergence (PSDTF)
}
\description{
The input data is assumed to be a three-dimensional array (M x M x N),
where each frontal slice is a positive semidefinite (PSD) matrix.
PSDTF approximates each PSD matrix as a weighted sum of K PSD basis matrices
using the LogDet divergence.
}
\usage{
PSDTF(X, rank=3, initH=NULL, initV=NULL,
  thr=1e-10, num.iter=100, verbose=FALSE)
}

\arguments{
\item{X}{
A three-dimensional array of size M x M x N, where each frontal slice
X[,,n] is a positive semidefinite matrix.
}
\item{rank}{
The number of PSD basis matrices (Default: 3).
}
\item{initH}{
The initial values of activation matrix H (K x N non-negative matrix)
(Default: NULL).
}
\item{initV}{
The initial values of PSD basis matrices, a list of K matrices each of
size M x M (Default: NULL).
}
\item{thr}{
The threshold for convergence based on relative change of LogDet divergence
(Default: 1e-10).
}
\item{num.iter}{
The number of maximum iterations (Default: 100).
}
\item{verbose}{
If verbose is TRUE, the iteration number and error are printed
(Default: FALSE).
}
}

\value{
  H : The activation coefficient matrix (K-rows and N-columns)

  V : A list of K PSD basis matrices, each of size M x M

  RecError : The LogDet divergence values at each iteration

  RelChange : The relative change of error at each iteration
}

\references{
K. Yoshii, R. Tomioka, D. Mochihashi, and M. Goto,
"Beyond NMF: Time-domain audio source separation without phase reconstruction",
ISMIR, pp.369-374, 2013.
}

\author{Koki Tsuyuzaki}

\examples{
if(interactive()){
  # Generate toy PSD tensor
  M <- 4
  N <- 10
  K <- 2
  V_true <- list()
  for(k in 1:K){
    A <- matrix(runif(M*M), M, M)
    V_true[[k]] <- A \%*\% t(A)
  }
  H_true <- matrix(runif(K*N), K, N)
  X <- array(0, dim=c(M, M, N))
  for(n in 1:N){
    for(k in 1:K){
      X[,,n] <- X[,,n] + H_true[k,n] * V_true[[k]]
    }
  }

  # Run PSDTF
  out <- PSDTF(X, rank=K, num.iter=20, verbose=TRUE)
}
}

\keyword{methods}
