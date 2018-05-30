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
matdata <- toyModel(model = "NMF")
tensordata1 <- toyModel(model = "CP")
tensordata2 <- toyModel(model = "Tucker")
}

\keyword{methods}