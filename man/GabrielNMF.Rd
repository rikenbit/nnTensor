\name{GabrielNMF}
\docType{methods}
\alias{GabrielNMF}

\title{
Gabriel-type Bi-Cross-Validation for Non-negative Matrix Factorization
}
\description{
The input data is assumed to be non-negative matrix.
GabrielNMF devides the input file into four matrices (A, B, C, and D)
and perform cross validation by the prediction of A from the matrices B, C, and D.
}
\usage{
GabrielNMF(X, J = 3, nx = 5, ny = 5,
 algorithm = c("Frobenius", "KL", "IS", "Pearson", "Hellinger", "Neyman",
    "Alpha", "Beta", "PGD", "HALS", "GCD"), Alpha = 1, Beta = 2,
    eta = 1e-04, thr1 = 1e-10, thr2 = 1e-10, tol = 1e-04,
    num.iter = 100, verbose = FALSE)
}

\arguments{
\item{X}{
The input matrix which has N-rows and M-columns.
}
\item{J}{
The number of low-dimension (J < {N, M}).
}
\item{nx}{
The number of hold-out in row-wise direction (2 < nx < N).
}
\item{ny}{
The number of hold-out in row-wise direction (2 < ny < M).
}
\item{algorithm}{
NMF algorithms. "Frobenius", "KL", "IS", "Pearson", "Hellinger", "Neyman",
"HALS", "PGD", "Alpha", "Beta", and "GCD" are available (Default: "KL").
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
When error change rate is lower than thr1, the iteration is terminated
(Default: 1E-10).
}
\item{thr2}{
If the minus-value is generated, replaced as thr2 (Default: 1E-10).
This value is used within the internal function .positive().
}
\item{tol}{
The tolerance parameter used in GCD algorithm.
}
\item{num.iter}{
The number of interation step (Default: 100).
}
\item{verbose}{
If verbose == TRUE, Error change rate is generated in console window.
}
}

\value{
  TestRecError : The reconstruction error calculated by Gabriel-style Bi-Cross
  Validation.
}

\references{
Art B. Owen et. al., (2009). Bi-Cross-Validation of the SVD and the Nonnegative Matrix Factorization. \emph{The Annals of Applied Statistics}
}
\author{Koki Tsuyuzaki}

\examples{
  if(interactive()){
    # Test data
    matdata <- toyModel(model = "NMF")

    # Bi-Cross-Validation
    BCV <- rep(0, length=5)
    names(BCV) <- 2:6
    for(j in seq(BCV)){
      print(j+1)
      BCV[j] <- mean(GabrielNMF(matdata, J=j+1, nx=2, ny=2)$TestRecError)
    }
    proper.rank <- as.numeric(names(BCV)[which(BCV == min(BCV))])
    
    # NMF
    out <- NMF(matdata, J=proper.rank)
  }
}
\keyword{methods}