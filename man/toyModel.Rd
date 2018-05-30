\name{toyModel}
\docType{methods}
\alias{toyModel}

\title{
Toy model data for using NMF, NTF, and NTD
}
\description{
The data is used for confirming the algorithm are properly working.
}
\usage{
toyModel(model = "CP")
}

\arguments{
  \item{model}{
Single character string is specified.
"NMF", "CP", and "Tucker" are available (Default: "CP").
}
}

\value{
If model is specified as "NMF", a matrix is generated.
Otherwise, a tensor is generated.
}
\author{Koki Tsuyuzaki}

\seealso{
\code{\link{NMF}}, \code{\link{NTF}}, \code{\link{NTD}}
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