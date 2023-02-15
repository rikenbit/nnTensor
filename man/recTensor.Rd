\name{recTensor}
\docType{methods}
\alias{recTensor}

\title{
Tensor Reconstruction from core tensor (S) and factor matrices (A)
}
\description{
Combined with plotTensor3D function and the result of NTF or NTD, the reconstructed tesor structure can be visullized.
}
\usage{
recTensor(S = NULL, A = NULL, idx = seq_along(dim(S)), reverse = FALSE)
}

\arguments{
  \item{S}{
K-order tensor object, which is defined as S4 class of rTensor package.
}
  \item{A}{
A list containing K factor matrices.
}
  \item{idx}{
The direction of mode-n muliplication (Default: 1:K).
For example idx=1 is defined. S x_1 A is calculated (x_1 : mode-1 multiplication).
 }
  \item{reverse}{
If reverse = TRUE, t(A[[n]]) is multiplicated to S (Default: FALSE).
  }
}

\value{
Tensor object, which is defined as S4 class of rTensor package.
}
\author{Koki Tsuyuzaki}

\seealso{
\code{\link[rTensor]{Tensor-class}}, \code{\link{NTF}}, \code{\link{NTD}}
}
\examples{
tensordata <- toyModel(model = "CP")
out <- NTF(tensordata, rank=3, algorithm="Beta-HALS", num.iter=2)
rec <- recTensor(out$S, out$A)
}

\keyword{methods}