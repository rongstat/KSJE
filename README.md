# Duo-Landmark Joint Embeddings for High-Dimensional Data


We develop a kernel spectral method that achieves joint embeddings of two independently observed high-dimensional noisy datasets. The proposed method automatically captures and leverages shared low-dimensional structures across datasets to enhance embedding quality. The obtained low-dimensional embeddings can be utilized for downstream tasks such as simultaneous clustering, data visualization, and denoising. The proposed method is justified by rigorous theoretical analysis, which guarantees its consistency in capturing the signal structures, and provides a geometric interpretation of the embeddings. 

The method is based on the paper:

Ding, X., and Ma, R. (2025+) Kernel Spectral Joint Embeddings for High-Dimensional Noisy Datasets Using Duo-Landmark Integral Operators. *Journal of the American Statistical Association* [https://arxiv.org/pdf/2203.00126.](https://arxiv.org/pdf/2405.12317)


# Content

The folder /code contains R scripts that reproduce our numerical simulations and analyses of the example single-cell omic datasets. We also provide the R function (see `main_function.R`) that implement our proposed algorithm.

`main_function.R`: script of R function that implements the proposed duo-landmark joint embedding algorithm.

`simulation_biclustering.R`: R script for simulation part 1: simultaneous clustering;

`simulation_manifold_learning.R`: R script for simulation part 2: noisy manifold learning;

`real_data_analysis.R`: R script for real data analysis;

`gamma1gamma2.R`: R script for embedding index set sensitivity analysis;

`unequal_sample_size.R`: R script for sample size unbalance sensitivity analysis;

`simulation_runtime.R`: R script for running time evaluations;

`simu_torus_sensitivity.R`: R script for bandwith hyperparameter and latent dimension sensitivity analyses;

`counter_expr.R`: R script for the counter examples presented in Section A of the Supplement.


# System Requirements

The method requires only a standard computer with enough RAM to support the operations defined by a user. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB
CPU: 4+ cores, 3.3+ GHz/core

The R implementation of the method is tested under R version 4.2.3, and requires the R packages: `RSpectra`,`Rfast`,`ggplot2`,`fossil`,`Seurat`,`SeuratData`,`uwot`,`HDF5Array`,`zellkonverter`,`cluster`.

# Quick Guide

For a **quick guide** to our method implementation in R, please check out [https://rongstat.github.io/KSJE_guide.io/tutorial.html](https://rongstat.github.io/KSJE_Guide.io/tutorial.html).

For further questions and inquiries, please contact Rong Ma (rongma@hsph.harvard.edu).
