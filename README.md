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
- Non-negative Matrix Factorization (NMF) : [Nonnegative Matrix and Tensor Factorizations, Andrzej CICHOCK, et. al., 2009](https://pdfs.semanticscholar.org/94cc/6daad548a03c6edb0351d686c2d4aa364634.pdf), [A Study on Efficient Algorithms for Nonnegative Matrix/Tensor Factorization, Keigo Kimura, 2017](https://eprints.lib.hokudai.ac.jp/dspace/bitstream/2115/65379/1/Keigo_Kimura.pdf)
- Simultaneous Non-negative Matrix Factorization (siNMF) : [Extracting Gene Expression Profiles Common to Colon and Pancreatic Adenocarcinoma using Simultaneous nonnegative matrix factorization, Liviu Badea, Pacific Symposium on Biocomputing, 13:279-290, 2009](https://psb.stanford.edu/psb-online/proceedings/psb08/badea.pdf), [Discovery of multi-dimensional modules by integrative analysis of cancer genomic data. Shihua Zhang, et al., Nucleic Acids Research, 40(19), 9379-9391, 2012](https://academic.oup.com/nar/article/40/19/9379/2414808), [Probabilistic Latent Tensor Factorization, International Conference on Latent Variable Analysis and Signal Separation, Y. Kenan Yilmaz et al., 346-353, 2010](http://papers.nips.cc/paper/3208-probabilistic-matrix-factorization.pdf)
- Joint Non-negative Matrix Factorization (jNMF) : [A non-negative matrix factorization method for detecting modules in heterogeneous omics multi-modal data, Zi Yang, et al., Bioinformatics, 32(1), 1-8, 2016](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btv544)
- Non-negative CP Decomposition (NTF)
   - α-Divergence (KL, Pearson, Hellinger, Neyman) / β-Divergence (KL, Frobenius, IS) : [Non-negative Tensor Factorization using Alpha and Beta Divergence, Andrzej CICHOCKI et. al., 2007](http://mlg.postech.ac.kr/static/publications/inter_conf/2007/icassp07_cichocki.pdf), [TensorKPD.R (gist of mathieubray)](https://gist.github.com/mathieubray/d83ce9c13fcb60f723f957c13ad85ac5)
   - Fast HALS : [Multi-way Nonnegative Tensor Factorization Using Fast Hierarchical Alternating Least Squares Algorithm (HALS), Anh Huy PHAN et. al., 2008](http://www.ieice.org/proceedings/NOLTA2008/articles/A1L-D3-Phan-2045.pdf)
   - α-HALS/β-HALS : [Fast Local Algorithms for Large Scale Nonnegative Matrix and Tensor Factorizations, Andrzej CICHOCKI et. al., 2008](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.214.6398&rep=rep1&type=pdf)
- Non-negative Tucker Decomposition (NTD)
   - KL, Frobenius : [Nonnegative Tucker Decomposition, Yong-Deok Kim et. al., 2007](https://pdfs.semanticscholar.org/f388/99be8ebd8b9aa7029b2b4f187dac4b04d816.pdf)
   - α-Divergence (KL, Pearson, Hellinger, Neyman) / β-Divergence (KL, Frobenius, IS) : [Nonneegative Tucker Decomposition With Alpha-Divergence, Yong-Deok Kim et. al., 2008](https://pdfs.semanticscholar.org/f01b/7354619f053863048217c58cc517def86aeb.pdf), [Fast and efficient algorithms for nonnegative Tucker decomposition, Anh Huy Phan, 2008](https://link.springer.com/chapter/10.1007/978-3-540-87734-9_88)
   - Fast HALS : [Extended HALS algorithm for nonnegative Tucker decomposition and its applications for multiway analysis and classification, Anh Hyu Phan et. al., 2011](https://www.sciencedirect.com/science/article/pii/S0925231211000427)
- Rank estimation of NMF
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

## License
Copyright (c) 2018 Koki Tsuyuzaki and Laboratory for Bioinformatics Research, RIKEN Center for Biosystems Dynamics Reseach
Released under the [Artistic License 2.0](http://www.perlfoundation.org/artistic_license_2_0).

## Authors
- Koki Tsuyuzaki
- Manabu Ishii
- Itoshi Nikaido