\name{plotTensor3D}
\docType{methods}
\alias{plotTensor3D}

\title{
Plot function for visualization of tensor data structure
}
\description{
Combined with recTensor function and the result of NTF or NTD, the reconstructed tesor structure can be visullized.
}
\usage{
plotTensor3D(X = NULL)
}

\arguments{
  \item{X}{
Tensor object, which is defined as S4 class of rTensor package.
  }
}
\author{Koki Tsuyuzaki}

\examples{
library("fields")

matdata <- toyModel(model = "NMF")
tensordata1 <- toyModel(model = "CP")
tensordata2 <- toyModel(model = "Tucker")

out1 <- NMF(matdata, J=5, algorithm="Frobenius")
out2 <- NTF(tensordata1, rank=4, algorithm="Beta-HALS")
out3 <- NTD(tensordata2, rank=c(5,5,5), algorithm="KL")

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