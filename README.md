# Duo-Landmark Joint Embeddings for High-Dimensional Data


We develop a kernel spectral method that achieves joint embeddings of two independently observed high-dimensional noisy datasets. The proposed method automatically captures and leverages shared low-dimensional structures across datasets to enhance embedding quality. The obtained low-dimensional embeddings can be utilized for downstream tasks such as simultaneous clustering, data visualization, and denoising. The proposed method is justified by rigorous theoretical analysis, which guarantees its consistency in capturing the signal structures, and provides a geometric interpretation of the embeddings. 

The method is based on the paper:

Ding, X., and Ma, R. (2025+) Kernel Spectral Joint Embeddings for High-Dimensional Noisy Datasets Using Duo-Landmark Integral Operators. *Journal of the American Statistical Association* [https://arxiv.org/pdf/2203.00126.](https://arxiv.org/pdf/2405.12317)


# Content

The folder /code contains R scripts that reproduce our numerical simulations and analyses of the example single-cell omic datasets. We also provide the R function (see `main_function.R`) that implement our proposed algorithm.

# System Requirements

The method requires only a standard computer with enough RAM to support the operations defined by a user. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB
CPU: 4+ cores, 3.3+ GHz/core

The R implementation of the method is tested under R version 4.2.3, and requires the R packages: `RSpectra`,`Rfast`,`ggplot2`,`fossil`,`Seurat`,`SeuratData`,`uwot`,`HDF5Array`,`zellkonverter`,`cluster`.
