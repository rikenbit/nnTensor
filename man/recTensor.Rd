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
recTensor(S = NULL, A = NULL, idx = 1:3, reverse = FALSE)
}

\arguments{
  \item{S}{
Tensor object, which is defined as S4 class of rTensor package.
}
  \item{A}{
A list containing three factor matrices.
}
  \item{idx}{
The direction of mode-n muliplication (Default: 1:3). For example idx=1 is defined. S x_1 A is calculated (x_1 : mode-1 multiplication).
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
library("fields")

matdata <- toyModel(model = "NMF")
tensordata1 <- toyModel(model = "CP")
tensordata2 <- toyModel(model = "Tucker")

out1 <- NMF(matdata, J=5, algorithm="Frobenius", num.iter=10)
out2 <- NTF(tensordata1, rank=4, algorithm="Beta-HALS", num.iter=10)
out3 <- NTD(tensordata2, rank=c(5,5,5), algorithm="KL", num.iter=10)

tmp <- tempdir()

png(filename=paste0(tmp, "/NMF.png"))
image.plot(matdata)
dev.off()

png(filename=paste0(tmp, "/NTF.png"))
plotTensor3D(recTensor(out2$S, out2$A))
dev.off()

png(filename=paste0(tmp, "/NTD.png"))
plotTensor3D(recTensor(out3$S, out3$A))
dev.off()
}

\keyword{methods}