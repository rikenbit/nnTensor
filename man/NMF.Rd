\name{NMF}
\docType{methods}
\alias{NMF}

\title{
Non-negative Matrix Factorization Algorithms (NMF)
}
\description{
The input data is assumed to be non-negative matrix.
NMF decompose the matrix to two low-dimensional factor matices.
This function is also used as initialization step of tensor decomposition (see also NTF and NTD).
}
\usage{
NMF(X, J = 3, algorithm = "KL", Alpha = 1,
  Beta = 2, eta = 1e-04, thr1 = 1e-10, thr2 = 1e-10,
  tol = 1e-04, num.iter = 100, viz = FALSE, figdir = NULL, verbose = FALSE)
}

\arguments{
  \item{X}{
The input Matrix which has N-rows and M-columns.
}
  \item{J}{
Number of low-dimension (J < N, M).
}
  \item{algorithm}{
NMF algorithms. "Frobenius", "KL", "IS", "Pearson", "Hellinger", "Neyman", "HALS", "PGD", "Alpha", "Beta", and "GCD" are available (Default: "KL").
}
  \item{Alpha}{
The parameter of Alpha-divergence.
}
  \item{Beta}{
The parameter of Beta-divergence.
}
  \item{eta}{
The stepsize for PGD algorithm (Default: 0.0001).
}
  \item{thr1}{
When error change rate is lower than thr1, the iteration is terminated (Default: 1E-10).
}
  \item{thr2}{
If the minus-value is generated, replaced as thr2 (Default: 1E-10). This value is used within the internal function .positive().
}
  \item{tol}{
The tolerance parameter used in GCD algorithm.
}
  \item{num.iter}{
The number of interation step (Default: 100).
}
  \item{viz}{
If viz == TRUE, internal reconstructed matrix can be visualized.
}
  \item{figdir}{
the directory for saving the figure, when viz == TRUE.
}
  \item{verbose}{
If verbose == TRUE, Error change rate is generated in console windos.
}
}

\value{
   U : A matrix which has N-rows and J-columns (J < N, M).
   V : A matrix which has M-rows and J-columns (J < N, M).
   RecError : The reconstruction error between data tensor and reconstructed tensor from U and V
   RelChange : The relative change of the error
}


\references{
Andrzej CICHOCK, et. al., (2009). Nonnegative Matrix and Tensor Factorizations. \emph{John Wiley & Sons, Ltd}

Keigo Kimura, (2017). A Study on Efficient Algorithms for Nonnegative Matrix/Tensor Factorization. \emph{Hokkaido University Collection of Scholarly and Academic Papers}
}
\author{Koki Tsuyuzaki}

\examples{
matdata <- toyModel(model = "NMF")
out <- NMF(matdata, J=2, num.iter=2)
}

\keyword{methods}