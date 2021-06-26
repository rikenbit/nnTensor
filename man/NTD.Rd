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
NTD(X, M=NULL, initS=NULL, initA=NULL, fixS=FALSE, fixA=FALSE,
  L1_A=1e-10, L2_A=1e-10,
  rank = c(3, 3, 3), modes = 1:3,
  algorithm = c("Frobenius", "KL", "IS", "Pearson", "Hellinger",
    "Neyman", "HALS", "Alpha", "Beta"), init = c("NMF", "ALS", "Random"),
  Alpha = 1,
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
  \item{initS}{
The initial values of core tensor which has J1, J2, and J3 dimensions
(Default: NULL).
}
  \item{initA}{
A list containing the initial values of multiple factor matrices
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
The number of low-dimension in each mode (J1, J2, J3, J1<I1, J2<I2, J3 < I3)
(Default: c(3,3,3)).
}
  \item{modes}{
The vector of the modes on whih to perform the decomposition
(Default: 1:3 <all modes>).
}
  \item{algorithm}{
NTD algorithms. "Frobenius", "KL", "IS", "Pearson", "Hellinger", "Neyman",
"HALS", "Alpha", and "Beta" are available (Default: "Frobenius").
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