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
kFoldMaskTensor(X, k=5, avoid.zero=TRUE, seeds=123)
}

\arguments{
  \item{X}{
An rTensor object.
  }
  \item{k}{
Number of split for k-fold cross validation.
  }
  \item{avoid.zero}{
If TRUE, only non-zero elements are splitted (Default: TRUE).
  }
  \item{seeds}{
Random seed to use for set.seed().
  }
}
\author{Koki Tsuyuzaki}

\examples{
tensordata <- toyModel(model = "CP")

Ms <- kFoldMaskTensor(tensordata, k=5, avoid.zero=TRUE, seeds=123)
}

\keyword{methods}