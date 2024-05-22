# Kernel Spectral Joint Embeddings for High-Dimensional Data


We develop a kernel spectral method that achieves joint nonlinear embeddings of two independently observed high-dimensional noisy datasets. The proposed method automatically captures and leverages the possibly shared low-dimensional structures across datasets to enhance embedding quality. The obtained low-dimensional embeddings can be utilized for many downstream tasks such as simultaneous clustering, data visualization, and denoising.  The proposed method is justified by rigorous theoretical analysis, which guarantees its statistical consistency in relation to the underlying signal structures, and provides  explicit geometric interpretations of the low-dimensional embeddings. 

The method is based on the paper:

Ding, X., and Ma, R. (2024+) Kernel Spectral Joint Embeddings for High-Dimensional Noisy Datasets Using Duo-Landmark Integral Operators. [https://arxiv.org/pdf/2203.00126.](https://arxiv.org/pdf/2405.12317)


# Content

The folder /data contains the example datasets analyzed in the paper.

The folder /code contains R scripts for numerical simulations and the analyses of the example datasets.

# System Requirements

The method requires only a standard computer with enough RAM to support the operations defined by a user. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB
CPU: 4+ cores, 3.3+ GHz/core

The R implementation of the method is tested under R version 4.2.3, and requires the R packages: `RSpectra`,`Rfast`,`ggplot2`,`fossil`,`Seurat`,`SeuratData`,`uwot`,`HDF5Array`,`zellkonverter`,`cluster`.
