\name{NTD}
\docType{methods}
\alias{NTD}

\title{
Non-negative Tucker Decomposition Algorithms (NTD)
}
\description{
The input data is assumed to be non-negative tensor.
NTD decompose the tensor to the dense core tensor (S) and low-dimensional
factor matices (A).
}
\usage{
NTD(X, M=NULL, pseudocount=.Machine$double.eps, initS=NULL, initA=NULL,
  fixS=FALSE, fixA=FALSE, L1_A=1e-10, L2_A=1e-10,
  rank = rep(3, length=length(dim(X))),
  modes = seq_along(dim(X)),
  algorithm = c("Frobenius", "KL", "IS", "Pearson", "Hellinger", "Neyman",
    "HALS", "Alpha", "Beta", "NMF"), init = c("NMF", "ALS", "Random"),
  nmf.algorithm = c("Frobenius", "KL", "IS", "Pearson", "Hellinger", "Neyman",
    "Alpha", "Beta", "ALS", "PGD", "HALS", "GCD", "Projected", "NHR", "DTPP",
    "Orthogonal", "OrthReg"),
  Alpha = 1,
  Beta = 2, thr = 1e-10, num.iter = 100, num.iter2 = 10, viz = FALSE,
  figdir = NULL, verbose = FALSE)
}

\arguments{
  \item{X}{
K-order input tensor which has I_1, I_2, ..., and I_K dimensions.
}
  \item{M}{
K-order mask tensor which has I_1, I_2, ..., and I_K dimensions. If the mask tensor has
missing values, specify the element as 0 (otherwise 1).
}
  \item{pseudocount}{
The pseudo count to avoid zero division, when the element is zero (Default: Machine Epsilon).
}
  \item{initS}{
The initial values of core tensor which has I_1, I_2, ..., and I_K dimensions
(Default: NULL).
}
  \item{initA}{
A list containing the initial values of K factor matrices
(A_k, <Ik*Jk>, k=1..K, Default: NULL).
}
  \item{fixS}{
Whether the core tensor S is updated in each iteration step (Default: FALSE).
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
The number of low-dimension in each mode (Default: 3 for each mode).
}
  \item{modes}{
The vector of the modes on which to perform the decomposition
(Default: 1:K <all modes>).
}
  \item{algorithm}{
NTD algorithms. "Frobenius", "KL", "IS", "Pearson", "Hellinger", "Neyman",
"HALS", "Alpha", "Beta", "NMF" are available (Default: "Frobenius").
}
  \item{nmf.algorithm}{
NMF algorithms, when the algorithm is "NMF".
"Frobenius", "KL", "IS", "Pearson", "Hellinger", "Neyman", "Alpha", "Beta",
"ALS", "PGD", "HALS", "GCD", "Projected", "NHR", "DTPP", "Orthogonal",
and "OrthReg" are available (Default: "Frobenius").
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
  \item{num.iter2}{
The number of NMF interation step, when the algorithm is "NMF" (Default: 10).
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
  S : K-order tensor object, which is defined as S4 class of rTensor package.
  A : A list containing K factor matrices.
  RecError : The reconstruction error between data tensor and reconstructed
  tensor from S and A.
  TrainRecError : The reconstruction error calculated by training set
  (observed values specified by M).
  TestRecError : The reconstruction error calculated by test set
  (missing values specified by M).
  RelChange : The relative change of the error.
}

\references{
Yong-Deok Kim et. al., (2007). Nonnegative Tucker Decomposition.
\emph{IEEE Conference on Computer Vision and Pattern Recognition}

Yong-Deok Kim et. al., (2008). Nonneegative Tucker Decomposition With
Alpha-Divergence. \emph{IEEE International Conference on Acoustics,
Speech and Signal Processing}

Anh Huy Phan, (2008). Fast and efficient algorithms for nonnegative
Tucker decomposition. \emph{Advances in Neural Networks - ISNN2008}

Anh Hyu Phan et. al. (2011). Extended HALS algorithm for nonnegative
Tucker decomposition and its applications for multiway analysis and classification.
\emph{Neurocomputing}
}
\author{Koki Tsuyuzaki}

\seealso{
\code{\link{plotTensor3D}}
}
\examples{
tensordata <- toyModel(model = "Tucker")
out <- NTD(tensordata, rank=c(2,2,2), algorithm="Frobenius",
  init="Random", num.iter=2)
}

\keyword{methods}