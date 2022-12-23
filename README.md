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
- **Non-negative Matrix Factorization (NMF)** : Nonnegative Matrix and Tensor Factorizations, Andrzej CICHOCK, et. al., 2009, A Study on Efficient Algorithms for Nonnegative Matrix/Tensor Factorization, Keigo Kimura, 2017
- **Projected NMF**
- **Nonnegative Hebbian Rule (NHR)**
- **Ding-Ti-Peg-Park (DTPP) algorithm**
- **(Column vector-wise) Orthogonal NMF**
  - Algorithms for Orthogonal Nonnegative Matrix Factorization, Seungjin Choi, 2008
- **(Column vector-wise) Orthogonality-regularized NMF**
  - Orthogonal matrix factorization enables integrative analysis of multiple RNA binding proteins, Martin Stražar, Marinka Žitnik, Blaž Zupan, Jernej Ule, Tomaž Curk, Bioinformatics, 15;32(10):1527-35, 2016
- **Non-negative Matrix Tri-Factorization (NMTF)** : Fast Optimization of Non-Negative Matrix Tri-Factorization: Supporting Information, Andrej Copar, et. al., PLOS ONE, 14(6), e0217994, 2019, Co-clustering by Block Value Decomposition, Bo Long et al., SIGKDD'05, 2005, Orthogonal Nonnegative Matrix Tri-Factorizations for Clustering, Chris Ding et. al., 12th ACM SIGKDD, 2006
- **Simultaneous Non-negative Matrix Factorization (siNMF)** : Extracting Gene Expression Profiles Common to Colon and Pancreatic Adenocarcinoma using Simultaneous nonnegative matrix factorization, Liviu Badea, Pacific Symposium on Biocomputing, 13:279-290, 2009, Discovery of multi-dimensional modules by integrative analysis of cancer genomic data. Shihua Zhang, et al., Nucleic Acids Research, 40(19), 9379-9391, 2012, Probabilistic Latent Tensor Factorization, International Conference on Latent Variable Analysis and Signal Separation, Y. Kenan Yilmaz et al., 346-353, 2010
- **Joint Non-negative Matrix Factorization (jNMF)** : A non-negative matrix factorization method for detecting modules in heterogeneous omics multi-modal data, Zi Yang, et al., Bioinformatics, 32(1), 1-8, 2016
- **Non-negative CP Decomposition (NTF)**
   - *α-Divergence (KL, Pearson, Hellinger, Neyman) / β-Divergence (KL, Frobenius, IS)* : Non-negative Tensor Factorization using Alpha and Beta Divergence, Andrzej CICHOCKI et. al., 2007, TensorKPD.R (gist of mathieubray)
   - *Fast HALS* : Multi-way Nonnegative Tensor Factorization Using Fast Hierarchical Alternating Least Squares Algorithm (HALS), Anh Huy PHAN et. al., 2008
   - *α-HALS/β-HALS* : Fast Local Algorithms for Large Scale Nonnegative Matrix and Tensor Factorizations, Andrzej CICHOCKI et. al., 2008
- **Non-negative Tucker Decomposition (NTD)**
   - *KL, Frobenius* : Nonnegative Tucker Decomposition, Yong-Deok Kim et. al., 2007
   - *α-Divergence (KL, Pearson, Hellinger, Neyman) / β-Divergence (KL, Frobenius, IS)* : Nonneegative Tucker Decomposition With Alpha-Divergence, Yong-Deok Kim et. al., 2008, Fast and efficient algorithms for nonnegative Tucker decomposition, Anh Huy Phan, 2008
   - *Fast HALS* : Extended HALS algorithm for nonnegative Tucker decomposition and its applications for multiway analysis and classification, Anh Hyu Phan et. al., 2011
- **Rank estimation of NMF**
	- Jean-Philippe Brunet. et. al., (2004). Metagenes and molecular pattern discovery using matrix factorization. PNAS
	- Xiaoxu Han. (2007). CANCER MOLECULAR PATTERN DISCOVERY BY SUBSPACE CONSENSUS KERNEL CLASSIFICATION
	- Attila Frigyesi. et. al., (2008). Non-Negative Matrix Factorization for the Analysis of Complex Gene Expression Data: Identification of Clinically Relevant Tumor Subtypes. Cancer Informatics
	- Haesun Park. et. al., (2019). Lecture 3: Nonnegative Matrix Factorization: Algorithms and Applications. SIAM Gene Golub Summer School, Aussois France, June 18, 2019
	- Chunxuan Shao. et. al., (2017). Robust classification of single-cell transcriptome data by nonnegative matrix factorization. Bioinformatics
	- Paul Fogel (2013). Permuted NMF: A Simple Algorithm Intended to Minimize the Volume of the Score Matrix
	- Philip M. Kim. et. al., (2003). Subsystem Identification Through Dimensionality Reduction of Large-Scale Gene Expression Data. Genome Research
	- Lucie N. Hutchins. et. al., (2008). Position-dependent motif characterization using non-negative matrix factorization. Bioinformatics
	- Patrik O. Hoyer (2004). Non-negative Matrix Factorization with Sparseness Constraints. Journal of Machine Learning 5
	- N. Fujita et al., (2018) Biomarker discovery by integrated joint non-negative matrix factorization and pathway signature analyses, Scientific Report
	- Art B. Owen et. al., (2009). Bi-Cross-Validation of the SVD and the Nonnegative Matrix Factorization. The Annals of Applied Statistics
- **Exponent term depending on Beta parameter**
  - M. Nakano et al., (2010). Convergence-guaranteed multiplicative algorithms for nonnegative matrix factorization with Beta-divergence. IEEE Workshop on Machine Learning for Signal Processing

## Authors
- Koki Tsuyuzaki
- Manabu Ishii
- Itoshi Nikaido