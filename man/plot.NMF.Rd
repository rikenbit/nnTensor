\name{plot.NMF}
\docType{methods}
\alias{plot.NMF}
\alias{plot}

\title{
Plot function of the result of NMF function
}
\description{
Only if J is specified as a vector longer than 1, this function will be active.
}
\usage{
plot(x, ...)
}

\arguments{
  \item{x}{
The result of NMF function (NMF class).
}
  \item{...}{
Optional parameter for plot.
}
}

\value{
A ggplot will be gerarated.
}

\references{
Jean-Philippe Brunet. et. al., (2004). Metagenes and molecular pattern discovery using matrix factorization. \emph{PNAS}

Xiaoxu Han. (2007). CANCER MOLECULAR PATTERN DISCOVERY BY SUBSPACE CONSENSUS KERNEL CLASSIFICATION

Attila Frigyesi. et. al., (2008). Non-Negative Matrix Factorization for the Analysis of Complex Gene Expression Data: Identification of Clinically Relevant Tumor Subtypes. \emph{Cancer Informatics}

Haesun Park. et. al., (2019). Lecture 3: Nonnegative Matrix Factorization: Algorithms and Applications. \emph{SIAM Gene Golub Summer School, Aussois France, June 18, 2019}

Chunxuan Shao. et. al., (2017). Robust classification of single-cell transcriptome data by nonnegative matrix factorization. \emph{Bioinformatics}

Paul Fogel (2013). Permuted NMF: A Simple Algorithm Intended to Minimize the Volume of the Score Matrix

Philip M. Kim. et. al., (2003). Subsystem Identification Through Dimensionality Reduction of Large-Scale Gene Expression Data. \emph{Genome Research}

Lucie N. Hutchins. et. al., (2008). Position-dependent motif characterization using non-negative matrix factorization. \emph{Bioinformatics}

Patrik O. Hoyer (2004). Non-negative Matrix Factorization with Sparseness Constraints. \emph{Journal of Machine Learning 5}
}
\author{Koki Tsuyuzaki}

\examples{
plot.NMF
}

\keyword{methods}