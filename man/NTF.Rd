\name{NTF}
\docType{methods}
\alias{NTF}

\title{
Non-negative CP Decomposition Algorithms (NTF)
}
\description{
The input data is assumed to be non-negative tensor.
NTF decompose the tensor to the diagonal core tensor (S) and low-dimensional
factor matices (A).
}
\usage{
NTF(X, M=NULL, initA=NULL, fixA=FALSE, L1_A=1e-10, L2_A=1e-10, rank = 3,
  algorithm = c("Frobenius", "KL", "IS", "Pearson", "Hellinger", "Neyman",
    "HALS", "Alpha-HALS", "Beta-HALS", "Alpha", "Beta"),
  init = c("NMF", "ALS", "Random"), Alpha = 1,
  Beta = 2, thr = 1e-10, num.iter = 100, viz = FALSE,
  figdir = NULL, verbose = FALSE)
}

\arguments{
  \item{X}{
The input tensor which has I1, I2, and I3 dimensions.
}
  \item{M}{
The mask tensor which has I1, I2, and I3 dimensions. If the mask tensor has
missing values, specify the element as 0 (otherwise 1).
}
  \item{initA}{
A list containing the initial values of multiple factor matrices
(A_k, <Ik*Jk>, k=1..K, Default: NULL).
}
  \item{fixA}{
Whether the factor matrices Ak are updated in each iteration step
(Default: FALSE).
}
  \item{L1_A}{
Paramter for L1 regularitation (Default: 1e-10).
This also works as small positive constant to prevent division by zero,
so should be set as 0.
}
  \item{L2_A}{
Paramter for L2 regularitation (Default: 1e-10).
}
  \item{rank}{
The number of low-dimension in each mode (J1=J2=J3, J1<I1, J2<I2, J3 < I3)
(Default: 3).
}
  \item{algorithm}{
NTF algorithms. "Frobenius", "KL", "IS", "Pearson", "Hellinger", "Neyman", "HALS", "Alpha-HALS", "Beta-HALS", "Alpha", and "Beta" are available
(Default: "Frobenius").
}
  \item{init}{
The initialization algorithms. "NMF", "ALS", and "Random" are available
(Default: "NMF").
}
  \item{Alpha}{
The parameter of Alpha-divergence.
}
  \item{Beta}{
The parameter of Beta-divergence.
}
  \item{thr}{
When error change rate is lower than thr1, the iteration is terminated
(Default: 1E-10).
}
  \item{num.iter}{
The number of interation step (Default: 100).
}
  \item{viz}{
If viz == TRUE, internal reconstructed tensor can be visualized.
}
  \item{figdir}{
the directory for saving the figure, when viz == TRUE (Default: NULL).
}
  \item{verbose}{
If verbose == TRUE, Error change rate is generated in console windos.
}
}

\value{
  S : Tensor object, which is defined as S4 class of rTensor package.
  A : A list containing three factor matrices.
  RecError : The reconstruction error between data tensor and reconstructed
  tensor from S and A.
  TrainRecError : The reconstruction error calculated by training set
  (observed values specified by M).
  TestRecError : The reconstruction error calculated by test set
  (missing values specified by M).
  RelChange : The relative change of the error.
}

\references{
Andrzej CICHOCKI et. al., (2007). Non-negative Tensor Factorization using
Alpha and Beta Divergence. \emph{IEEE ICASSP 2007}

Anh Huy PHAN et. al., (2008). Multi-way Nonnegative Tensor Factorization
Using Fast Hierarchical Alternating Least Squares Algorithm (HALS). \emph{NOLTA2008}

Andrzej CICHOCKI et. al., (2008). Fast Local Algorithms for Large Scale
Nonnegative Matrix and Tensor Factorizations.
\emph{IEICE Transactions on Fundamentals of Electronics, Communications
and Computer Sciences}
}
\author{Koki Tsuyuzaki}

\seealso{
\code{\link{plotTensor3D}}
}
\examples{
tensordata <- toyModel(model = "CP")
out <- NTF(tensordata, rank=3, algorithm="Beta-HALS", num.iter=2)
}

\keyword{methods}