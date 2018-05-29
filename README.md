# nnTensor
R package for Non-negative Tensor Decomposition

Installation
======
~~~~
git clone https://github.com/rikenbit/nnTensor/
cd nnTensor
R CMD INSTALL nnTensor
~~~~

References
======
- Non-negative Matrix Factorization (NMF) : [Nonnegative Matrix and Tensor Factorizations, Andrzej CICHOCK, et. al., 2009](https://pdfs.semanticscholar.org/94cc/6daad548a03c6edb0351d686c2d4aa364634.pdf), [A Study on Efficient Algorithms for Nonnegative Matrix/Tensor Factorization, Keigo Kimura, 2017](https://eprints.lib.hokudai.ac.jp/dspace/bitstream/2115/65379/1/Keigo_Kimura.pdf)
- Non-negative CP Decomposition (NTF)
   - α-Divergence (KL, Pearson, Hellinger, Neyman) / β-Divergence (KL, Frobenius, IS) : [Non-negative Tensor Factorization using Alpha and Beta Divergence, Andrzej CICHOCKI et. al., 2007](http://mlg.postech.ac.kr/static/publications/inter_conf/2007/icassp07_cichocki.pdf), [TensorKPD.R (gist of mathieubray)](https://gist.github.com/mathieubray/d83ce9c13fcb60f723f957c13ad85ac5)
   - Fast HALS : [Multi-way Nonnegative Tensor Factorization Using Fast Hierarchical Alternating Least Squares Algorithm (HALS), Anh Huy PHAN et. al., 2008](http://www.ieice.org/proceedings/NOLTA2008/articles/A1L-D3-Phan-2045.pdf)
   - α-HALS/β-HALS : [Fast Local Algorithms for Large Scale Nonnegative Matrix and Tensor Factorizations, Andrzej CICHOCKI et. al., 2008](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.214.6398&rep=rep1&type=pdf)
- Non-negative Tucker Decomposition (NTD)
   - KL, Frobenius : [Nonnegative Tucker Decomposition, Yong-Deok Kim et. al., 2007](https://pdfs.semanticscholar.org/f388/99be8ebd8b9aa7029b2b4f187dac4b04d816.pdf)
   - α-Divergence (KL, Pearson, Hellinger, Neyman) / β-Divergence (KL, Frobenius, IS) : [Nonneegative Tucker Decomposition With Alpha-Divergence, Yong-Deok Kim et. al., 2008](https://pdfs.semanticscholar.org/f01b/7354619f053863048217c58cc517def86aeb.pdf), [Fast and efficient algorithms for nonnegative Tucker decomposition, Anh Huy Phan, 2008](https://link.springer.com/chapter/10.1007/978-3-540-87734-9_88)
   - Fast HALS : [Extended HALS algorithm for nonnegative Tucker decomposition and its applications for multiway analysis and classification Author links open overlay panel, Anh Hyu Phan et. al., 2011](https://www.sciencedirect.com/science/article/pii/S0925231211000427)

## License
Copyright (c) 2018 Koki Tsuyuzaki and Laboratory for Bioinformatics Research, RIKEN Center for Biosystems Dynamics Reseach
Released under the [Artistic License 2.0](http://www.perlfoundation.org/artistic_license_2_0).

## Authors
- Koki Tsuyuzaki
- Itoshi Nikaido
