# Discovery of spatial patterns with Hidden Markov random field

Use HMRF to dissect spatial patterns in smFISH data


## Introduction

`smfishHmrf` package provides tools for spatial pattern discovery
using multivariate Gaussian models and hidden Markov random field fitted by
expectation-maximization. We focus on seqFISH data.

[single-molecule (sm)FISH](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3126657/)
 is used to spatially profile cells at single-cell resolution. 
[seqFISH](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3418883/), 
a type of smFISH utilizing sequential multiplexed barcoding schemes, 
can profile hundreds of genes' transcript count for individual cells with high sensitivity.

Using gene expression from seqFISH, we can infer the spatial pattern that might
exist in the spatial profile. This is analogous to image segmentation in computer
vision. In the single-cell imaging, spatial pattern looks like a contiguous
domain of cells, sharing expression of genes, and may suggest something about 
the cell's local environment. This spatial state of the cell is unknown and 
needs to be identified from seqFISH data.

As the spatial states are assumed to be independent, we can use a mixture of multivariate
Gaussian distribution to model the observed expression with parameters estimated by the EM
algorithm. Since nearby cells tend to be of the same state, a **Markov
random field model** (in this case **Pott's model**) can be used to capture
the spatial similarity of cells by making homogeneous the relationship
between neighboring cells.

Typical package (such as [ref1](http://www.tandfonline.com/doi/abs/10.1198/jasa.2011.ap09529)) 
runs HMRF on magnetic resonance imaging (MRI) data and is restricted to 1-channel (grayscale) 
and fixed to 4-neighbor regular grid of pixels.
In our case, the package `smfishHmrf` is more general:

1. It extends simple 1D Gaussian to
multivariate Gaussian with multidimensional mean and covariance matrices estimated by
EM. 

2. As the cells do not fall on a regular grid, we adapt to a local neighborhood
graph instead of a grid. This allows for a more flexible specification of number of neighbors
and neighbor structure. 

To illustrate the method, we have tested it on
mouse brain visual cortex seqFISH imaged cells.


## Installation

See [INSTALL.md](INSTALL.md)
