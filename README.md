[![DOI](https://zenodo.org/badge/135140554.svg)](https://zenodo.org/badge/latestdoi/135140554)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/nnTensor)](
https://cran.r-project.org/package=nnTensor)
[![Downloads](https://cranlogs.r-pkg.org/badges/nnTensor)](https://CRAN.R-project.org/package=nnTensor)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/nnTensor?color=orange)](https://CRAN.R-project.org/package=nnTensor)
[![:name status badge](https://rikenbit.r-universe.dev/badges/:name)](https://rikenbit.r-universe.dev)
[![:registry status badge](https://rikenbit.r-universe.dev/badges/:registry)](https://rikenbit.r-universe.dev)
[![:total status badge](https://rikenbit.r-universe.dev/badges/:total)](https://rikenbit.r-universe.dev)
[![nnTensor status badge](https://rikenbit.r-universe.dev/badges/nnTensor)](https://rikenbit.r-universe.dev)
![GitHub Actions](https://github.com/rikenbit/nnTensor/actions/workflows/build_test_push.yml/badge.svg)
[![status](https://joss.theoj.org/papers/b8cc3029f784ee95c45831467b3b2f74/status.svg)](https://joss.theoj.org/papers/b8cc3029f784ee95c45831467b3b2f74)

# nnTensor
R package for Non-negative Tensor Decomposition

Installation
======
~~~~
git clone https://github.com/rikenbit/nnTensor/
R CMD INSTALL nnTensor
~~~~
or type the code below in the R console window
~~~~
library(devtools)
devtools::install_github("rikenbit/nnTensor")
~~~~

References
======
- **Non-negative Matrix Factorization (NMF)**
  - Lee, D. and Seung, H. Learning the parts of objects by non-negative matrix factorization. Nature 401, 788–791 (1999)
  - Cichocki, A. et al., Nonnegative Matrix and Tensor Factorizations, Wiley, 2009
  - Kimura, K. A Study on Efficient Algorithms for Nonnegative Matrix/Tensor Factorization, Ph.D. Thesis, 2017
- **Projective NMF/Nonnegative Hebbian Rule (NHR)/Ding-Ti-Peg-Park (DTPP) algorithm/(Column vector-wise) Orthogonal NMF**
  - Choi, S. Algorithms for Orthogonal Nonnegative Matrix Factorization, IEEE World Congress on Computational Intelligence, 1828-1832, 2008
- **(Column vector-wise) Orthogonality-regularized NMF**
  - Stražar, M. et al., Orthogonal matrix factorization enables integrative analysis of multiple RNA binding proteins, Bioinformatics, 32(10), 1527-35, 2016
- **Non-negative Matrix Tri-Factorization (NMTF)**
  - Copar, A. et al., Fast Optimization of Non-Negative Matrix Tri-Factorization: Supporting Information, PLOS ONE, 14(6), e0217994, 2019
  - Long, B. et al., Co-clustering by Block Value Decomposition, SIGKDD'05, 635–640, 2005
  - Ding, C. et al., Orthogonal Nonnegative Matrix Tri-Factorizations for Clustering, 12th ACM SIGKDD'06, 126–135, 2006
- **Simultaneous Non-negative Matrix Factorization (siNMF)**
  - Badea, L. Extracting Gene Expression Profiles Common to Colon and Pancreatic Adenocarcinoma using Simultaneous nonnegative matrix factorization, Pacific Symposium on Biocomputing, 279-290, 2008
  - Zhang, S. et al., Discovery of multi-dimensional modules by integrative analysis of cancer genomic data. Nucleic Acids Research, 40(19), 9379-9391, 2012
  - Yilmaz, Y. K. et al., Probabilistic Latent Tensor Factorization, IVA/ICA 2010, 346-353, 2010
- **Joint Non-negative Matrix Factorization (jNMF)**
  - Zi, Yang, et al., A non-negative matrix factorization method for detecting modules in heterogeneous omics multi-modal data, Bioinformatics, 32(1), 1-8, 2016
- **Non-negative CP Decomposition (NTF)**
   - *α-Divergence (KL, Pearson, Hellinger, Neyman) / β-Divergence (KL, Frobenius, IS)*
     - Cichocki, A. et al., Non-negative Tensor Factorization using Alpha and Beta Divergence, ICASSP '07, III-1393-III-1396, 2007
     - [mathieubray/TensorKPD.R](https://gist.github.com/mathieubray/d83ce9c13fcb60f723f957c13ad85ac5)
   - *Fast HALS*
     - Phan, A. H. et al.,  Multi-way Nonnegative Tensor Factorization Using Fast Hierarchical Alternating Least Squares Algorithm (HALS), NOLTA 2008, 2008
   - *α-HALS/β-HALS*
     - Cichocki, A. et al., Fast Local Algorithms for Large Scale Nonnegative Matrix and Tensor Factorizations, IEICE Transactions, 92-A, 708-721, 2009
- **Non-negative Tucker Decomposition (NTD)**
   - *Frobenius/KL*
     - Kim, Y.-D. et al., Nonnegative Tucker Decomposition, IEEE CVPR, 1-8, 2007
   - *α-Divergence (KL, Pearson, Hellinger, Neyman) / β-Divergence (KL, Frobenius, IS)*
     - Kim, Y.-D. et al., Nonneegative Tucker Decomposition with Alpha-Divergence, 2008
     - Phan, A. H. et al., Fast and efficient algorithms for nonnegative Tucker decomposition, ISNN 2008, 772-782, 2008
   - *Fast HALS*
     - Phan, A. H. et al., Extended HALS algorithm for nonnegative Tucker decomposition and its applications for multiway analysis and classification, Neurocomputing, 74(11), 1956-1969, 2011
- **Rank estimation of NMF**
	- Brunet, J.-P. et al., Metagenes and molecular pattern discovery using matrix factorization. PNAS, 101(12), 4164-4169, 2004
	- Han, X. Cancer Molecular Pattern Discovery by Subspace Consensus Kernel Classification. CSB 2007, 6, 55-65, 2007
	- Frigyesi, A. et al., Non-Negative Matrix Factorization for the Analysis of Complex Gene Expression Data: Identification of Clinically Relevant Tumor Subtypes. Cancer Informatics, 2008
	- Park, H. et al., Lecture 3: Nonnegative Matrix Factorization: Algorithms and Applications. SIAM Gene Golub Summer School, 2019
	- Shao, C. et al., Robust classification of single-cell transcriptome data by nonnegative matrix factorization. Bioinformatics, 33(2), 235-242, 2017
	- Fogel, P., Permuted NMF: A Simple Algorithm Intended to Minimize the Volume of the Score Matrix, arXiv, 2013
	- Kim, P. M. et al., Subsystem Identification Through Dimensionality Reduction of Large-Scale Gene Expression Data. Genome Research, 13(7), 1706-1718, 2003
	- Hutchins, L. N. et al., Position-dependent motif characterization using non-negative matrix factorization. Bioinformatics, 24(23), 2684-2690, 2008
	- Hoyer, P. O., Non-negative Matrix Factorization with Sparseness Constraints. JMLR 5, 1457-1469, 2004
	- Fujita, N. et al., Biomarker discovery by integrated joint non-negative matrix factorization and pathway signature analyses, Scientific Report, 8(1), 9743, 2018
	- Owen, A. B. et al., Bi-Cross-Validation of the SVD and the Nonnegative Matrix Factorization. The Annals of Applied Statistics, 3(2), 564-594, 2009
- **Exponent term depending on Beta parameter**
  - Nakano, M. et al., Convergence-guaranteed multiplicative algorithms for nonnegative matrix factorization with Beta-divergence. IEEE MLSP, 283-288, 2010

## Authors
- Koki Tsuyuzaki
- Manabu Ishii
- Itoshi Nikaido