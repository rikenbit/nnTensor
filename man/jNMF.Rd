\name{jNMF}
\docType{methods}
\alias{jNMF}

\title{
Joint Non-negative Matrix Factorization Algorithms (jNMF)
}
\description{
The input data objects are assumed to be non-negative matrices.
jNMF decompose the matrices to two low-dimensional factor matices simultaneously.
}
\usage{
jNMF(X, J = 3, w=NULL, algorithm = "KL", p=1,
    thr = 1e-10, num.iter = 100,
    viz = FALSE, figdir = NULL, verbose = FALSE)
}

\arguments{
  \item{X}{
A list containing input matrices (X_k, <N*Mk>, k=1..K).
}
  \item{J}{
Number of low-dimension (J < N, Mk).
}
  \item{w}{
Weight vector (Default: NULL)
}
  \item{algorithm}{
Divergence between X and X_bar. "Frobenius", "KL", and "IS" are available (Default: "KL").
}
  \item{p}{
The parameter of Probabilistic Latent Tensor Factorization (p=0: Frobenius, p=1: KL, p=2: IS)
}
  \item{thr}{
When error change rate is lower than thr, the iteration is terminated (Default: 1E-10).
}
  \item{num.iter}{
The number of interation step (Default: 100).
}
  \item{viz}{
If viz == TRUE, internal reconstructed matrix can be visualized.
}
  \item{figdir}{
the directory for saving the figure, when viz == TRUE.
}
  \item{verbose}{
If verbose == TRUE, Error change rate is generated in console windos.
}
}

\value{
   W : A matrix which has N-rows and J-columns (J < N, Mk).
   V : A list which has multiple elements containing N-rows and J-columns (J < N, Mk).
   H : A list which has multiple elements containing Mk-rows and J-columns matrix (J < N, Mk).
   RecError : The reconstruction error between data matrix and reconstructed matrix from W and H
   RelChange : The relative change of the error
}

\references{
Liviu Badea, (2008) Extracting Gene Expression Profiles Common to Colon and Pancreatic Adenocarcinoma using Simultaneous nonnegative matrix factorization. \emph{Pacific Symposium on Biocomputing} 13:279-290

Shihua Zhang, et al. (2012) Discovery of multi-dimensional modules by integrative analysis of cancer genomic data. \emph{Nucleic Acids Research} 40(19), 9379-9391

Zi Yang, et al. (2016) A non-negative matrix factorization method for detecting modules in heterogeneous omics multi-modal data, \emph{Bioinformatics} 32(1), 1-8

Y. Kenan Yilmaz et al., (2010) Probabilistic Latent Tensor Factorization, \emph{International Conference on Latent Variable Analysis and Signal Separation} 346-353
}
\author{Koki Tsuyuzaki}

\examples{
matdata <- toyModel(model = "siNMF_Hard")
out <- jNMF(matdata, J=2, num.iter=2)
}

\keyword{methods}