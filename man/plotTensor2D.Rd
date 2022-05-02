\name{plotTensor2D}
\docType{methods}
\alias{plotTensor2D}

\title{
Plot function for visualization of matrix data structure
}
\description{
Combined with recTensor function and the result of NTF or NTD, the reconstructed tensor structure can be visullized.
}
\usage{
plotTensor2D(X = NULL, method=c("sd", "mad"),
    sign=c("positive", "negative", "both"), thr=2)
}

\arguments{
  \item{X}{
Matrix object.
  }
  \item{method}{
Cutoff method to focus on large/small value in the tensor data (Default: "sd").
  }
  \item{sign}{
Direction to cutoff the large/small value in the tensor data (Default: "positive").
  }
  \item{thr}{
Threshold of cutoff method (Default: 2).
  }
}
\author{Koki Tsuyuzaki}

\examples{
tensordata <- toyModel(model = "CP")

out <- NTF(tensordata, rank=3, num.iter=2)

tmp <- tempdir()

png(filename=paste0(tmp, "/NTF.png"))
plotTensor2D(out$A[[1]])
dev.off()
}

\keyword{methods}