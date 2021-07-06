\name{siNMF}
\docType{methods}
\alias{siNMF}

\title{
Simultaneous Non-negative Matrix Factorization Algorithms (siNMF)
}
\description{
The input data objects are assumed to be non-negative matrices.
siNMF decompose the matrices to two low-dimensional factor
matices simultaneously.
}
\usage{
siNMF(X, M=NULL, pseudocount=1e-10, initW=NULL, initH=NULL, fixW=FALSE, fixH=FALSE,
  L1_W=1e-10, L1_H=1e-10, L2_W=1e-10, L2_H=1e-10, J = 3,
  w=NULL, algorithm = c("Frobenius", "KL", "IS", "PLTF"), p=1,
    thr = 1e-10, num.iter = 100,
    viz = FALSE, figdir = NULL, verbose = FALSE)
}

\arguments{
  \item{X}{
A list containing the input matrices (X_k, <N*Mk>, k=1..K).
}
  \item{M}{
A list containing the mask matrices (X_k, <N*Mk>, k=1..K). If the input matrix
has missing values, specify the element as 0 (otherwise 1).
}
  \item{pseudocount}{
The pseudo count to avoid zero division, when the element is zero (Default: 1e-10).
}
  \item{initW}{
The initial values of factor matrix W, which has N-rows and J-columns
(Default: NULL).
}
  \item{initH}{
A list containing the initial values of multiple factor matrices
(H_k, <Mk*J>, k=1..K, Default: NULL).
}
  \item{fixW}{
Whether the factor matrix W is updated in each iteration step (Default: FALSE).
}
  \item{fixH}{
Whether the factor matrices Hk are updated in each iteration step
(Default: FALSE).
}
  \item{L1_W}{
Paramter for L1 regularitation (Default: 1e-10).
This also works as small positive constant to prevent division by zero,
so should be set as 0.
}
  \item{L1_H}{
Paramter for L1 regularitation (Default: 1e-10).
This also works as small positive constant to prevent division by zero,
so should be set as 0.
}
  \item{L2_W}{
Paramter for L2 regularitation (Default: 1e-10).
}
  \item{L2_H}{
Paramter for L2 regularitation (Default: 1e-10).
}
  \item{J}{
Number of low-dimension (J < N, Mk).
}
  \item{w}{
Weight vector (Default: NULL)
}
  \item{algorithm}{
Divergence between X and X_bar. "Frobenius", "KL", and "IS" are available
(Default: "KL").
}
  \item{p}{
The parameter of Probabilistic Latent Tensor Factorization
(p=0: Frobenius, p=1: KL, p=2: IS)
}
  \item{thr}{
When error change rate is lower than thr, the iteration is terminated
(Default: 1E-10).
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
  H : A list which has multiple elements containing Mk-rows and
  J-columns matrix (J < N, Mk).
  RecError : The reconstruction error between data matrix and
  reconstructed matrix from W and H.
  TrainRecError : The reconstruction error calculated by training set
  (observed values specified by M).
  TestRecError : The reconstruction error calculated by test set
  (missing values specified by M).
  RelChange : The relative change of the error.
}

\references{
Liviu Badea, (2008) Extracting Gene Expression Profiles Common to Colon
and Pancreatic Adenocarcinoma using Simultaneous nonnegative matrix
factorization. \emph{Pacific Symposium on Biocomputing} 13:279-290

Shihua Zhang, et al. (2012) Discovery of multi-dimensional modules by
integrative analysis of cancer genomic data. \emph{Nucleic Acids Research}
40(19), 9379-9391

Zi Yang, et al. (2016) A non-negative matrix factorization method for
detecting modules in heterogeneous omics multi-modal data,
\emph{Bioinformatics} 32(1), 1-8

Y. Kenan Yilmaz et al., (2010) Probabilistic Latent Tensor Factorization,
\emph{International Conference on Latent Variable Analysis and Signal Separation}
346-353

N. Fujita et al., (2018) Biomarker discovery by integrated joint non-negative matrix factorization and pathway signature analyses, \emph{Scientific Report}
}
\author{Koki Tsuyuzaki}

\examples{
matdata <- toyModel(model = "siNMF_Easy")
out <- siNMF(matdata, J=2, num.iter=2)
}

\keyword{methods}