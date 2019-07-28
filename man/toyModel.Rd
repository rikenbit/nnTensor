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
toyModel(model = "CP", seeds=123)
}

\arguments{
  \item{model}{
Single character string is specified.
"NMF", "CP", and "Tucker" are available (Default: "CP").
}
  \item{seeds}{
Random number for setting set.seeds in the function (Default: 123).
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
matdata <- toyModel(model = "NMF", seeds=123)
tensordata1 <- toyModel(model = "CP", seeds=123)
tensordata2 <- toyModel(model = "Tucker", seeds=123)
}

\keyword{methods}