\name{NMTF}
\docType{methods}
\alias{NMTF}

\title{
Non-negative Matrix Tri-Factorization Algorithms (NMTF)
}
\description{
The input data is assumed to be non-negative matrix.
NMTF decompose the matrix to three low-dimensional factor matices.
}
\usage{
NMTF(X, pseudocount=.Machine$double.eps,
    initU=NULL, initS=NULL, initV=NULL,
    fixU=FALSE, fixS=FALSE, fixV=FALSE,
    L1_U=1e-10, L1_S=1e-10, L1_V=1e-10,
    L2_U=1e-10, L2_S=1e-10, L2_V=1e-10,
    orthU=FALSE, orthV=FALSE,
    rank = c(3, 4),
    algorithm = c("Frobenius", "KL", "IS", "ALS", "PG", "COD", "Beta"),
    Beta = 2, root = FALSE, thr = 1e-10, num.iter = 100,
    viz = FALSE, figdir = NULL, verbose = FALSE)
}

\arguments{
\item{X}{
The input matrix which has N-rows and M-columns.
}
\item{pseudocount}{
The pseudo count to avoid zero division, when the element is zero (Default: Machine Epsilon).
}
\item{initU}{
The initial values of factor matrix U, which has N-rows and J1-columns
(Default: NULL).
}
\item{initS}{
The initial values of factor matrix S, which has J1-rows and J2-columns
(Default: NULL).
}
\item{initV}{ The initial values of factor matrix V, which
has M-rows and J2-columns (Default: NULL).
}
\item{fixU}{ Whether the factor matrix U is updated in each iteration step
(Default: FALSE).
}
\item{fixS}{ Whether the factor matrix S is updated in each iteration step
(Default: FALSE).
}
\item{fixV}{
Whether the factor matrix V is updated in each iteration step
(Default: FALSE).
}
  \item{L1_U}{
Paramter for L1 regularitation (Default: 1e-10).
}
\item{L1_S}{
Paramter for L1 regularitation (Default: 1e-10).
}
\item{L1_V}{
Paramter for L1 regularitation (Default: 1e-10).
}
\item{L2_U}{
Paramter for L2 regularitation (Default: 1e-10).
}
\item{L2_S}{
Paramter for L2 regularitation (Default: 1e-10).
}
\item{L2_V}{
Paramter for L2 regularitation (Default: 1e-10).
}
\item{orthU}{
Whether the column vectors of matrix U are orthogonalized (Default: FALSE).
}
\item{orthV}{
Whether the column vectors of matrix V are orthogonalized (Default: FALSE).
}
\item{rank}{
The number of low-dimension (J1 (< N) and J2 (< M)) (Default: c(3,4)).
}
\item{algorithm}{
NMTF algorithms. "Frobenius", "KL", "IS", "ALS", "PG", "COD", and "Beta" are available (Default: "Frobenius").
}
\item{Beta}{
The parameter of Beta-divergence (Default: 2, which means "Frobenius").
}
\item{root}{
Whether square root is calculed in each iteration (Default: FALSE).
}
\item{thr}{
When error change rate is lower than thr, the iteration is terminated (Default: 1E-10).
}
\item{num.iter}{
The number of interation step (Default: 100).
}
\item{viz}{
If viz == TRUE, internal reconstructed matrix can be visualized.
}
\item{figdir}{
The directory for saving the figure, when viz == TRUE.
}
\item{verbose}{
If verbose == TRUE, Error change rate is generated in console window.
}
}

\value{
  U : A matrix which has N-rows and J1-columns (J1 < N).
  S : A matrix which has J1-rows and J2-columns.
  V : A matrix which has M-rows and J2-columns (J2 < M).
  RecError : The reconstruction error between data tensor and reconstructed
  tensor from U and V.
  RelChange : The relative change of the error.
  algorithm: algorithm specified.
}

\references{
Fast Optimization of Non-Negative Matrix Tri-Factorization: Supporting Information, Andrej Copar, et. al., \emph{PLOS ONE}, 14(6), e0217994, 2019

Co-clustering by Block Value Decomposition, Bo Long et al., \emph{SIGKDD'05}, 2005

Orthogonal Nonnegative Matrix Tri-Factorizations for Clustering, Chris Ding et. al., \emph{12th ACM SIGKDD}, 2006
}
\author{Koki Tsuyuzaki}

\examples{
  if(interactive()){
    # Test data
    matdata <- toyModel(model = "NMF")

    # Simple usage
    out <- NMTF(matdata, rank=c(4,4))
  }
}
\keyword{methods}