\name{kFoldMaskTensor}
\docType{methods}
\alias{kFoldMaskTensor}

\title{
Mask tensors generator to perform k-fold cross validation
}
\description{
The output multiple mask tensors can be immediately specified as the argument M for NTF() or NTD().
}
\usage{
kFoldMaskTensor(X, k=3, seeds=123, sym=FALSE)
}

\arguments{
  \item{X}{
An rTensor object.
  }
  \item{k}{
Number of split for k-fold cross validation (Default: 3).
  }
  \item{seeds}{
Random seed to use for set.seed() (Default: 123).
  }
  \item{sym}{
Data will be dropped symmetrically (available only when matrix is specified, Default: FALSE).
  }
}
\author{Koki Tsuyuzaki}

\examples{
tensordata <- toyModel(model = "CP")

str(kFoldMaskTensor(tensordata, k=5))
}

\keyword{methods}