\name{NMF}
\docType{methods}
\alias{NMF}

\title{
Non-negative Matrix Factorization Algorithms (NMF)
}
\description{
The input data is assumed to be non-negative matrix.
NMF decompose the matrix to two low-dimensional factor matices.
This function is also used as initialization step of tensor decomposition
(see also NTF and NTD).
}
\usage{
NMF(X, M=NULL, pseudocount=.Machine$double.eps, initU=NULL, initV=NULL,
  fixU=FALSE, fixV=FALSE,
  L1_U=1e-10, L1_V=1e-10, L2_U=1e-10, L2_V=1e-10, J = 3,
  rank.method=c("all", "ccc", "dispersion", "rss", "evar", "residuals",
    "sparseness.basis", "sparseness.coef", "sparseness2.basis",
    "sparseness2.coef", "norm.info.gain.basis", "norm.info.gain.coef",
    "singular",  "volume", "condition"), runtime=30,
  algorithm = c("Frobenius", "KL", "IS", "Pearson", "Hellinger", "Neyman",
    "Alpha", "Beta", "ALS", "PGD", "HALS", "GCD", "Projected", "NHR", "DTPP",
    "Orthogonal", "OrthReg", "t", "TV-GNMF"),
  Alpha = 1, Beta = 2, nu = 2,
  eta = 1e-04, thr1 = 1e-10, thr2 = 1e-10, tol = 1e-04,
  num.iter = 100,
  L_graph_U = NULL, L_graph_V = NULL,
  lambda_graph_U = 0, lambda_graph_V = 0,
  graph.method = c("Vicus", "LEM", "HLLE"),
  graph.K = 10, graph.alpha = 0.9,
  lambda_TV = 0, beta_TV = 0,
  viz = FALSE, figdir = NULL, verbose = FALSE)
}

\arguments{
\item{X}{
The input matrix which has N-rows and M-columns.
}
\item{M}{
The mask matrix which has N-rows and M-columns. If the input matrix has
missing values, specify the elements as 0 (otherwise 1).
}
  \item{pseudocount}{
The pseudo count to avoid zero division, when the element is zero (Default: Machine Epsilon).
}
\item{initU}{
The initial values of factor matrix U, which has N-rows and J-columns
(Default: NULL).
}
\item{initV}{ The initial values of factor matrix V, which
has M-rows and J-columns (Default: NULL).
}
\item{fixU}{ Whether the factor matrix U is updated in each iteration step
(Default: FALSE).
}
\item{fixV}{
Whether the factor matrix V is updated in each iteration step
(Default: FALSE).
}
  \item{L1_U}{
Paramter for L1 regularitation (Default: 1e-10).
This also works as small positive constant to prevent division by zero,
so should be set as 0.
}
  \item{L1_V}{
Paramter for L1 regularitation (Default: 1e-10).
This also works as small positive constant to prevent division by zero,
so should be set as 0.
}
  \item{L2_U}{
Paramter for L2 regularitation (Default: 1e-10).
}
  \item{L2_V}{
Paramter for L2 regularitation (Default: 1e-10).
}
\item{J}{
The number of low-dimension (J < \{N, M\}). If a numerical vector is specified
(e.g. 2:6), the appropriate rank is estimated.
}
\item{rank.method}{
The rank estimation method (Default: "all"). Only if the J option is
specified as a numerical vector longer than two, this option will be active.
}
\item{runtime}{
The number of trials to estimate rank (Default: 10).
}
\item{algorithm}{
NMF algorithms. "Frobenius", "KL", "IS", "Pearson", "Hellinger", "Neyman",
"Alpha", "Beta", "ALS", "PGD", "HALS", "GCD", "Projected", "NHR", "DTPP",
"Orthogonal", "OrthReg", "t", and "TV-GNMF" are available
(Default: "Frobenius").
}
\item{Alpha}{
The parameter of Alpha-divergence.
}
\item{Beta}{
The parameter of Beta-divergence.
}
\item{nu}{
The degrees of freedom parameter for t-NMF algorithm (Default: 2).
Larger values approach Gaussian NMF. Only used when algorithm = "t".
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
\item{L_graph_U}{
A pre-computed graph Laplacian matrix for rows (N x N).
If NULL and lambda_graph_U > 0, it is computed automatically using
Vicus::graphMatrix (Default: NULL).
}
\item{L_graph_V}{
A pre-computed graph Laplacian matrix for columns (M x M).
If NULL and lambda_graph_V > 0, it is computed automatically using
Vicus::graphMatrix (Default: NULL).
}
\item{lambda_graph_U}{
The weight for graph regularization on U (Default: 0).
When set to a positive value, adds graph Laplacian regularization to the
denominator of the multiplicative update rule for U.
}
\item{lambda_graph_V}{
The weight for graph regularization on V (Default: 0).
When set to a positive value, adds graph Laplacian regularization to the
denominator of the multiplicative update rule for V.
}
\item{graph.method}{
The method for computing the graph matrix via Vicus::graphMatrix.
"Vicus", "LEM", and "HLLE" are available (Default: "Vicus").
}
\item{graph.K}{
The number of nearest neighbors for graph construction (Default: 10).
}
\item{graph.alpha}{
The alpha parameter for Vicus::graphMatrix (Default: 0.9).
}
\item{lambda_TV}{
The weight for Total Variation regularization on V.
Only used when algorithm = "TV-GNMF" (Default: 0).
}
\item{beta_TV}{
The weight for graph regularization in TV-GNMF algorithm.
Only used when algorithm = "TV-GNMF" (Default: 0).
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
  U : A matrix which has N-rows and J-columns (J < \{N, M\}).
  V : A matrix which has M-rows and J-columns (J < \{N, M\}).
  J : The number of dimension (J < \{N, M\}).
  RecError : The reconstruction error between data tensor and reconstructed
  tensor from U and V.
  TrainRecError : The reconstruction error calculated by training set
  (observed values specified by M).
  TestRecError : The reconstruction error calculated by test set
  (missing values specified by M).
  RelChange : The relative change of the error.
  Trial : All the results of the trials to estimate the rank.
  Runtime : The number of the trials to estimate the rank.
  RankMethod : The rank estimation method.
}

\references{
Andrzej CICHOCK, et. al., (2009). Nonnegative Matrix and Tensor Factorizations.
\emph{John Wiley & Sons, Ltd}

Keigo Kimura, (2017). A Study on Efficient Algorithms for Nonnegative Matrix/
Tensor Factorization.
\emph{Hokkaido University Collection of Scholarly and Academic Papers}

K. Yoshii, K. Itoyama, and M. Goto (2016). Student's t Nonnegative Matrix
Factorization and Positive Semidefinite Tensor Factorization for
Single-Channel Audio Source Separation. \emph{ICASSP}, pp.51-55.

C. Leng, Z. Chen, G. Cai, I. Cheng, Z. Xiong, J. Tian, and A. Basu (2016).
Total Variation Constrained Graph Regularized NMF for Medical Image
Registration. \emph{IEEE ICIP}.

D. Cai, X. He, J. Han, and T. S. Huang (2011). Graph Regularized Nonnegative
Matrix Factorization for Data Representation.
\emph{IEEE Transactions on PAMI}, 33(8), pp.1548-1560.
}
\author{Koki Tsuyuzaki}

\examples{
  if(interactive()){
    # Test data
    matdata <- toyModel(model = "NMF")

    # Simple usage
    out <- NMF(matdata, J=5)

    # Rank estimation mode (single method)
    out2 <- NMF(matdata, J=2:10, rank.method="ccc", runtime=3)
    plot(out2)

    # Rank estimation mode (all method)
    out3 <- NMF(matdata, J=2:10, rank.method="all", runtime=10)
    plot(out3)

    # t-NMF (robust to outliers)
    out4 <- NMF(matdata, J=5, algorithm="t", nu=2)

    # Graph-regularized NMF
    out5 <- NMF(matdata, J=5, algorithm="Frobenius",
      lambda_graph_V=0.1, graph.K=5)

    # TV-GNMF (Total Variation + Graph regularization)
    out6 <- NMF(matdata, J=5, algorithm="TV-GNMF",
      lambda_TV=0.1, beta_TV=0.1,
      lambda_graph_V=0.1, graph.K=5)
  }
}
\keyword{methods}
