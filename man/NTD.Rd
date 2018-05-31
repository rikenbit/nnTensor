\name{NTD}
\docType{methods}
\alias{NTD}

\title{
Non-negative Tucker Decomposition Algorithms (NTD)
}
\description{
The input data is assumed to be non-negative tensor.
NTD decompose the tensor to the dense core tensor (S) and low-dimensional factor matices (A).
}
\usage{
NTD(X, rank = c(3, 3, 3), algorithm = "KL", init = "NMF", Alpha = 1,
  Beta = 2, thr = 1e-10, num.iter = 100, viz = FALSE,
  figdir = ".", verbose = FALSE)
}

\arguments{
  \item{X}{
The input tensor which has I1, I2, and I3 dimensions.
}
  \item{rank}{
The number of low-dimension in each mode (J1, J2, J3, J1<I1, J2<I2, J3 < I3) (Default: c(3,3,3)).
}
  \item{algorithm}{
NTD algorithms. "Frobenius", "KL", "IS", "Pearson", "Hellinger", "Neyman", "HALS", "Alpha", and "Beta" are available (Default: "Frobenius").
}
  \item{init}{
The initialization algorithms. "NMF", "ALS", and "Random" are available (Default: "NMF").
}
  \item{Alpha}{
The parameter of Alpha-divergence.
}
  \item{Beta}{
The parameter of Beta-divergence.
}
  \item{thr}{
When error change rate is lower than thr1, the iteration is terminated (Default: 1E-10).
}
  \item{num.iter}{
The number of interation step (Default: 100).
}
  \item{viz}{
If viz == TRUE, internal reconstructed tensor can be visualized.
}
  \item{figdir}{
the directory for saving the figure, when viz == TRUE.
}
  \item{verbose}{
If verbose == TRUE, Error change rate is generated in console windos.
}
}

\value{
  S : Tensor object, which is defined as S4 class of rTensor package.
  A : A list containing three factor matrices.
}

\references{
Yong-Deok Kim et. al., (2007). Nonnegative Tucker Decomposition. \emph{IEEE Conference on Computer Vision and Pattern Recognition}

Yong-Deok Kim et. al., (2008). Nonneegative Tucker Decomposition With Alpha-Divergence. \emph{IEEE International Conference on Acoustics, Speech and Signal Processing}

Anh Huy Phan, (2008). Fast and efficient algorithms for nonnegative Tucker decomposition. \emph{Advances in Neural Networks - ISNN2008}

Anh Hyu Phan et. al. (2011). Extended HALS algorithm for nonnegative Tucker decomposition and its applications for multiway analysis and classification. \emph{Neurocomputing}
}
\author{Koki Tsuyuzaki}

\seealso{
\code{\link{plotTensor3D}}
}
\examples{
tensordata <- toyModel(model = "Tucker")
out <- NTD(tensordata, rank=c(2,2,2), algorithm="Frobenius", init="Random", num.iter=2)
}

\keyword{methods}